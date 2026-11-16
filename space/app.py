import importlib.util
import os
from pathlib import Path

import gradio as gr
import uvicorn
from fastapi import FastAPI
from fastapi.responses import RedirectResponse
from fastapi.staticfiles import StaticFiles

here = Path(__file__).resolve().parent
root = here if (here / "app" / "app.py").exists() else here.parent

spec = importlib.util.spec_from_file_location("brca_tool", root / "app" / "app.py")
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)
tool.explorer_url = "/explorer/"

server = FastAPI()


@server.get("/explorer", include_in_schema=False)
def explorer_slash():
    return RedirectResponse("/explorer/")


@server.get("/tool", include_in_schema=False)
def tool_slash():
    return RedirectResponse("/tool/")


@server.get("/landing", include_in_schema=False)
@server.get("/landing/", include_in_schema=False)
def landing_old():
    return RedirectResponse("/")


server.mount("/explorer", StaticFiles(directory=root / "explorer", html=True), name="explorer")
server = gr.mount_gradio_app(
    server,
    tool.demo,
    path="/tool",
    theme=tool.theme,
    css=tool.css_text,
    head=tool.font_link + tool.head_script,
    footer_links=[],
)
server.mount("/", StaticFiles(directory=root / "landing", html=True), name="landing")

if __name__ == "__main__":
    uvicorn.run(server, host="0.0.0.0", port=int(os.environ.get("PORT", "7860")))
