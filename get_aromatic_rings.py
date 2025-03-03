from amsr.encode import FromSmiles
import requests
from bs4 import BeautifulSoup
import json


def get_names(url="https://en.wikipedia.org/wiki/Simple_aromatic_ring"):
    response = requests.get(url)
    soup = BeautifulSoup(response.content, "html.parser")

    # Find the table with a caption containing "Table of simple aromatic rings"
    table = None
    for candidate in soup.find_all("table"):
        caption = candidate.find("caption")
        if caption and "Table of simple aromatic rings" in caption.get_text():
            table = candidate
            break

    names = set()
    if table:
        rows = table.find_all("tr")
        for row in rows:
            # Look only at rows that contain data cells
            cells = row.find_all("td")
            if not cells:
                continue
            for cell in cells:
                # In each cell, there are usually two <a> tags:
                # one linking to the image (File:) and one linking to the actual article.
                # We skip the image link by checking that the href doesn't start with "/wiki/File:"
                links = cell.find_all("a")
                for link in links:
                    href = link.get("href", "")
                    if href.startswith("/wiki/File:"):
                        continue
                    title = link.get("title")
                    if title:
                        names.add(title)
                        break  # Use the first valid name per cell
    return names


def get_smiles(name):
    url = f"https://en.wikipedia.org/wiki/{name}"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, "html.parser")

    # Find the anchor whose text is "SMILES" (case-insensitive)
    anchor = soup.find("a", string=lambda t: t and "smiles" in t.lower())
    if anchor:
        # Look for the next <ul> element with class "mw-collapsible-content"
        collapsible = anchor.find_next("ul", class_="mw-collapsible-content")
        if collapsible:
            # Assume the SMILES string is within the first <li><div> inside the collapsible content.
            li = collapsible.find("li")
            if li:
                inner_div = li.find("div")
                if inner_div:
                    return inner_div.get_text(strip=True)
    return None


smi = {n.lower(): get_smiles(n) for n in sorted(get_names())}

tab = {}
for name, s in smi.items():
    if s is None:
        continue
    n = 0
    for a in sorted(
        {FromSmiles(s, randomize=True, useGroups=False) for i in range(100000)}
    ):
        k = "(@)" if name == "benzene" and n == 0 else f"({name}_{n})"
        tab[k] = [a]
        n += 1

with open("aromatic_rings.json", "w") as f:
    json.dump(tab, f, indent=2)
    f.write("\n")
