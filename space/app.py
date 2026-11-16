import importlib.util
from pathlib import Path

from fastapi.responses import RedirectResponse
from fastapi.staticfiles import StaticFiles
from starlette.routing import Mount, Route

here = Path(__file__).resolve().parent
root = here if (here / "app" / "app.py").exists() else here.parent

try:
    import spaces

    @spaces.GPU
    def unused():
        return None
except ImportError:
    spaces = None

spec = importlib.util.spec_from_file_location("brca_tool", root / "app" / "app.py")
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)
tool.explorer_url = "/explorer/"


def explorer_slash(request):
    return RedirectResponse("/explorer/")


def attach(app):
    app.router.routes.insert(0, Route("/explorer", explorer_slash))
    app.router.routes.insert(1, Mount("/explorer", app=StaticFiles(directory=root / "explorer", html=True), name="explorer"))


if __name__ == "__main__":
    app, local_url, share_url = tool.demo.launch(
        theme=tool.theme,
        css=tool.css_text,
        head=tool.font_link + tool.head_script,
        footer_links=[],
        prevent_thread_lock=True,
    )
    attach(app)
    tool.demo.block_thread()
