use thiserror::Error;

#[derive(Debug, Error)]
pub enum SemError {
    #[error("syntax error: {0}")]
    SyntaxError(String),

    #[error("model did not converge after {iterations} iterations")]
    ConvergenceFailed { iterations: usize },

    #[error("singular weight matrix")]
    SingularWeight,

    #[error("model is not identified: {0}")]
    NotIdentified(String),

    #[error("matrix error: {0}")]
    MatrixError(String),
}
