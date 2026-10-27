import queue
import threading
from tkinter import messagebox


def run_in_background(widget, work, on_success, on_progress=None, on_error=None):
    """
    Run `work(progress)` in a background thread and deliver its result to
    `on_success(result)` in the GUI thread. Tkinter is not thread-safe, so the
    worker never touches widgets; it calls `progress(message)` instead, which
    is forwarded to `on_progress` in the GUI thread.
    """
    events = queue.Queue()

    def target():
        try:
            result = work(lambda message: events.put(("progress", message)))
            events.put(("done", result))
        except Exception as e:
            events.put(("error", e))

    def poll():
        while True:
            try:
                kind, value = events.get_nowait()
            except queue.Empty:
                widget.after(100, poll)
                return

            if kind == "progress":
                if on_progress:
                    on_progress(value)
            elif kind == "done":
                on_success(value)
                return
            else:
                (on_error or _show_error)(value)
                return

    threading.Thread(target=target, daemon=True).start()
    widget.after(100, poll)


def _show_error(error):
    messagebox.showerror("Error", str(error))
