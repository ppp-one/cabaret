from dataclasses import dataclass


@dataclass
class Focuser:
    """
    A simple focuser model for simulating the effect of defocus on image quality.

    Parameters
    ----------
    position : float
        The current focus position of the focuser.
    best_position : float
        The optimal focus position for best image quality.
    scale : float
        The scale factor that determines how quickly defocus degrades seeing.
    max_seeing_multiplier : float
        The maximum factor by which seeing can be increased due to defocus.

    Examples
    --------
    >>> from cabaret.focuser import Focuser
    >>> focuser = Focuser()
    >>> focuser.seeing_multiplier
    1.0
    >>> focuser.position = 10100
    >>> focuser.seeing_multiplier
    2.0

    Notes
    -----
    The `seeing_multiplier` property quantifies how much the atmospheric seeing is
    increased due to the focuser being away from the best focus position. The further
    the position is from `best_position`, the larger the multiplier, up to
    `max_seeing_multiplier`.
    """

    position: float = 10_000
    best_position: float = 10_000
    scale: float = 100
    max_seeing_multiplier: float = 5.0

    @property
    def seeing_multiplier(
        self,
    ) -> float:
        """Factor by which the seeing is increased due to defocus."""
        offset = abs(self.position - self.best_position)
        return min(1 + offset / self.scale, self.max_seeing_multiplier)
