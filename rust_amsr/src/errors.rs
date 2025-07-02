use thiserror::Error;

pub type AMSRResult<T> = Result<T, AMSRError>;

#[derive(Error, Debug)]
pub enum AMSRError {
    #[error("Invalid SMILES: {0}")]
    InvalidSmiles(String),
    
    #[error("Invalid AMSR: {0}")]
    InvalidAMSR(String),
    
    #[error("Invalid atom symbol: {0}")]
    InvalidAtomSymbol(String),
    
    #[error("Invalid bond: {0}")]
    InvalidBond(String),
    
    #[error("Parse error: {0}")]
    ParseError(String),
    
    #[error("Graph error: {0}")]
    GraphError(String),
    
    #[error("Valence error: {0}")]
    ValenceError(String),
    
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
} 
