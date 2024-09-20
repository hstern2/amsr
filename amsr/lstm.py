import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
from torch.nn.utils.rnn import pad_sequence

# Define special tokens
PAD_TOK = 0
EOS_TOK = 1


def _device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class LSTMModel(nn.Module):
    def __init__(self, d_model=128, nhid=256, nlayers=2, dropout=0.3):
        super(LSTMModel, self).__init__()
        self.d_model = d_model
        self.nhid = nhid
        self.nlayers = nlayers
        self.dropout = dropout

    def _init_model(self):
        self.embedding = nn.Embedding(self.vocab_size, self.d_model)
        self.lstm = nn.LSTM(
            input_size=self.d_model,
            hidden_size=self.nhid,
            num_layers=self.nlayers,
            dropout=self.dropout,
            batch_first=True,
        )
        self.fc_out = nn.Linear(self.nhid, self.vocab_size)
        self.to(_device())

    def forward(self, src):
        embedded = self.embedding(src)
        lstm_out, _ = self.lstm(embedded)
        output = self.fc_out(lstm_out)
        return output

    def _get_vocab(self, seqs):
        self.vocab_size = 2
        self.index_for_token = {}
        self.token_for_index = {}
        for s in seqs:
            for t in s:
                if t not in self.index_for_token:
                    self.index_for_token[t] = self.vocab_size
                    self.token_for_index[self.vocab_size] = t
                    self.vocab_size += 1

    def train_and_save_model(
        self,
        seqs,
        model_path,
        num_epochs=200,
        patience=10,
        learning_rate=0.001,
        weight_decay=1e-5,
        batch_size=64,
        validation_split=0.2,
    ):
        self._get_vocab(seqs)
        self._init_model()
        seqs_as_tt = [
            torch.tensor(
                [self.index_for_token[t] for t in s] + [EOS_TOK], device=_device()
            )
            for s in seqs
        ]
        padded_sequences = pad_sequence(
            seqs_as_tt, batch_first=True, padding_value=PAD_TOK
        )
        targets = padded_sequences[:, 1:]
        pad = torch.full(
            (targets.shape[0], 1), PAD_TOK, dtype=targets.dtype, device=_device()
        )
        padded_targets = torch.cat((targets, pad), dim=1)

        dataset = TensorDataset(padded_sequences, padded_targets)
        train_size = int((1 - validation_split) * len(dataset))
        val_size = len(dataset) - train_size
        train_dataset, val_dataset = random_split(dataset, [train_size, val_size])

        train_data_loader = DataLoader(
            train_dataset, batch_size=batch_size, shuffle=True
        )
        val_data_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

        self.train()
        criterion = nn.CrossEntropyLoss(ignore_index=PAD_TOK)
        optimizer = optim.Adam(
            self.parameters(), lr=learning_rate, weight_decay=weight_decay
        )

        best_val_loss = float("inf")
        epochs_no_improve = 0

        for epoch in range(num_epochs):
            self.train()  # training phase
            total_train_loss = 0
            for src, tgt in train_data_loader:
                src, tgt = src.to(_device()), tgt.to(_device())
                optimizer.zero_grad()
                output = self(src)
                loss = criterion(output.view(-1, self.vocab_size), tgt.view(-1))
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.parameters(), max_norm=1.0)
                optimizer.step()
                total_train_loss += loss.item()
            avg_train_loss = total_train_loss / len(train_data_loader)

            self.eval()  # validation phase
            total_val_loss = 0
            with torch.no_grad():
                for src, tgt in val_data_loader:
                    src, tgt = src.to(_device()), tgt.to(_device())
                    output = self(src)
                    loss = criterion(output.view(-1, self.vocab_size), tgt.view(-1))
                    total_val_loss += loss.item()
            avg_val_loss = total_val_loss / len(val_data_loader)
            print(
                f"Epoch {epoch}, Train Loss: {avg_train_loss}, Validation Loss: {avg_val_loss}",
                flush=True,
            )

            if avg_val_loss < best_val_loss:
                best_val_loss = avg_val_loss
                epochs_no_improve = 0
            else:
                epochs_no_improve += 1

            if epochs_no_improve >= patience:
                print(f"Early stopping at epoch {epoch}", flush=True)
                break

        torch.save(
            {
                "token_for_index": self.token_for_index,
                "index_for_token": self.index_for_token,
                "model_state_dict": self.state_dict(),
                "d_model": self.d_model,
                "nhid": self.nhid,
                "nlayers": self.nlayers,
                "vocab_size": self.vocab_size,
            },
            model_path,
        )
        print("Model and configuration saved to", model_path, flush=True)

    @classmethod
    def from_saved_model(cls, model_path):
        checkpoint = torch.load(model_path, map_location=_device(), weights_only=True)
        model = cls(
            d_model=checkpoint["d_model"],
            nhid=checkpoint["nhid"],
            nlayers=checkpoint["nlayers"],
        )
        model.vocab_size = checkpoint["vocab_size"]
        model.token_for_index = checkpoint["token_for_index"]
        model.index_for_token = checkpoint["index_for_token"]
        model._init_model()
        model.load_state_dict(checkpoint["model_state_dict"])
        model.eval()
        model.to(_device())
        return model

    def generate_tokens(self, start_input, max_length=20, temperature=1.0):

        assert temperature > 0

        self.eval()
        generated_sequence = torch.tensor(
            [self.index_for_token.get(t, PAD_TOK) for t in start_input],
            device=_device(),
        ).unsqueeze(0)

        for _ in range(max_length - 1):
            with torch.no_grad():
                output = self(generated_sequence)
                logits = (
                    output[:, -1, :] / temperature
                )  # Scale the logits by the temperature
                probabilities = torch.softmax(logits, dim=-1)
                next_token = torch.multinomial(
                    probabilities, num_samples=1
                )  # Sample from the probability distribution
                t = next_token.item()
                if t == PAD_TOK:
                    continue
                elif t == EOS_TOK:
                    break
                generated_sequence = torch.cat((generated_sequence, next_token), dim=1)

        n = len(start_input)
        return [
            start_input[_] if _ < n else self.token_for_index[i]
            for _, i in enumerate(generated_sequence.squeeze(0).tolist())
        ]

    def generate(self, start_input, max_length=20, temperature=1.0):
        return "".join(self.generate_tokens(start_input, max_length, temperature))

    def replace_token(self, seq, position, temperature=1.0):

        assert temperature > 0
        assert position >= 0 and position < len(seq)

        self.eval()
        seq_tensor = torch.tensor(
            [self.index_for_token.get(t, PAD_TOK) for t in seq], device=_device()
        ).unsqueeze(0)

        with torch.no_grad():
            output = self(seq_tensor)
            logits = output[:, position, :] / temperature
            probabilities = torch.softmax(logits, dim=-1)
            new_token_index = torch.multinomial(probabilities, num_samples=1).item()
            while new_token_index in (PAD_TOK, EOS_TOK):
                new_token_index = torch.multinomial(probabilities, num_samples=1).item()
            return (
                seq[:position]
                + [self.token_for_index[new_token_index]]
                + seq[position + 1 :]
            )


if __name__ == "__main__":
    seqs = [["A", "B", "C"], ["B", "C", "D", "E"], ["C", "A", "D"]]
    pth = "tmp_model.pth"
    model = LSTMModel()
    model.train_and_save_model(seqs, pth)
    loaded_model = LSTMModel.from_saved_model(pth)
    print(loaded_model.generate(["A", "B"]))
    print(loaded_model.replace_token(["B", "D", "D", "E"], 1))
