use thiserror::Error;

/// Errors returned by the [`MajorTomGrid`](crate::MajorTomGrid) API.
#[derive(Debug, Error)]
pub enum GridError {
    /// The grid spacing `d` was invalid (must be strictly greater than zero).
    #[error("grid spacing must be positive")]
    InvalidSpacing,

    /// The supplied cell ID was malformed (e.g. shorter than the geohash
    /// precision of 11 characters).
    #[error("cell ID must be at least 11 characters")]
    InvalidCellId,

    /// No cell matching the supplied ID could be reconstructed on the grid.
    #[error("no cell found with ID {0}")]
    CellNotFound(String),
}
