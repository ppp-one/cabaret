import concurrent.futures
import socket


def has_internet(host="duckduckgo.com", port=53, timeout=3):
    """Check for internet connectivity by attempting to connect to a reliable host."""
    try:
        socket.setdefaulttimeout(timeout)
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.connect((host, port))
        return True
    except OSError:
        return False


def has_tap_source(tap_source, timeout: int = 5) -> bool:
    """Check if a TAP service successfully processes a minimal query within `timeout` s.

    Probing /availability is not sufficient — the Gaia Archive returns HTTP 200
    from that endpoint even when queries time out or return 500 errors.  Running
    a real (tiny) query is the only reliable signal.
    """
    from astroquery.utils.tap.core import TapPlus

    from cabaret.queries import _TAP_CONFIG

    cfg = _TAP_CONFIG[tap_source]
    query = f"SELECT TOP 1 {cfg['ra']} FROM {cfg['from']}"

    def _probe():
        tap = TapPlus(url=tap_source.value)
        job = tap.launch_job(query)
        job.get_results()
        return True

    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            return executor.submit(_probe).result(timeout=timeout)
    except Exception:
        return False
