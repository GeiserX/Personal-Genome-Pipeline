"""A plugin for tests/test_pgp_live.py: one route, one run file glob and one state hook."""


def register(app):
    app.add_argument("--fixture-word", default="hello")

    def hello(req, m):
        app.log("fixture", m.group(1), status=200)
        req.send_json({"word": app.config.fixture_word, "run": m.group(1), "dir": bool(req.run_dir(m.group(1)))})

    app.route("GET", r"/fixture/([A-Za-z0-9._-]+)", hello)
    app.run_file("extra/*.txt")
    app.state_hook("fixture", lambda run_dir, state, events: {"events": len(events), "steps": len(state["steps"])})
