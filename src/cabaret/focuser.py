from dataclasses import dataclass


@dataclass
class Focuser:
    """
    A simple focuser model for simulating the effect of defocus on image quality.
    """

    position: float = 10_000
    """The current focus position of the focuser."""

    best_position: float = 10_000
    """The optimal focus position for best image quality."""

    scale: float = 100
    """The scale factor that determines how quickly defocus degrades seeing."""

    max_seeing_multiplier: float = 5.0
    """The maximum factor by which seeing can be increased due to defocus."""

    @property
    def seeing_multiplier(
        self,
    ) -> float:
        """Factor by which the seeing is increased due to defocus."""
        offset = abs(self.position - self.best_position)
        return min(1 + offset / self.scale, self.max_seeing_multiplier)
