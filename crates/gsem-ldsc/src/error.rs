use thiserror::Error;

#[derive(Debug, Error)]
pub enum LdscError {
    #[error("insufficient SNPs: need at least {min}, got {got}")]
    InsufficientSnps { min: usize, got: usize },

    #[error("LDSC regression failed to converge")]
    ConvergenceFailed,

    #[error("singular matrix in LDSC regression")]
    SingularMatrix,

    #[error("invalid input: {0}")]
    InvalidInput(String),
}
