use thiserror::Error;

#[derive(Debug, Clone, Error)]
pub enum MatrixError {
    #[error("matrix is not square: {rows}x{cols}")]
    NotSquare { rows: usize, cols: usize },

    #[error("matrix is not symmetric")]
    NotSymmetric,

    #[error("dimension mismatch: expected {expected}, got {got}")]
    DimensionMismatch { expected: usize, got: usize },

    #[error("singular matrix")]
    SingularMatrix,

    #[error("nearest PD did not converge after {iterations} iterations")]
    NearPdNotConverged { iterations: usize },
}
