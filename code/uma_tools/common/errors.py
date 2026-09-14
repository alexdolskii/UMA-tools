"""Input validation errors shared by the sequential UMA workflows."""


class ValidationError(ValueError):
    """Carry an input error and optional diagnostic records."""

    def __init__(self, message: str, details: list | None = None) -> None:
        super().__init__(message)
        self.details = details or []
