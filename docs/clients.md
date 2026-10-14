# Connect an MCP client

proteomics-mcp runs locally as an MCP server over stdio. Every client below starts the same command and only registers it differently:

```text
uvx proteomics-mcp@0.1.9 serve
```

`uvx` comes with [uv](https://docs.astral.sh/uv/getting-started/installation/). It installs the pinned release from PyPI into a cached environment on first use, so the first start takes a minute. Pin the version, as shown here, so that results are reproducible. On macOS, install the OpenMP runtime once for local MS2PIP prediction: `brew install libomp`.

Do not run `proteomics-mcp serve` in a terminal by itself; it waits for JSON-RPC messages from a client.

## Claude Desktop

Recommended: download `proteomics-mcp-v0.1.9.mcpb` from the [GitHub release](https://github.com/UadKLab/proteomics-mcp/releases/tag/v0.1.9) and open it with Claude Desktop. The extension manages its own Python environment. In its settings, choose an output directory; input files must be in folders Claude Desktop can read.

Alternatively, without the extension, add the server to `claude_desktop_config.json` (*Settings → Developer → Edit Config*) and restart Claude Desktop:

```json
{
  "mcpServers": {
    "proteomics-mcp": {
      "command": "uvx",
      "args": ["proteomics-mcp@0.1.9", "serve"]
    }
  }
}
```

## Claude Code

```bash
claude mcp add proteomics-mcp -- uvx proteomics-mcp@0.1.9 serve
```

Add `--scope project` to share the server with a repository through its `.mcp.json`. Check it with `claude mcp list`.

## Codex CLI

```bash
codex mcp add proteomics-mcp -- uvx proteomics-mcp@0.1.9 serve
```

or add it to `~/.codex/config.toml`:

```toml
[mcp_servers.proteomics-mcp]
command = "uvx"
args = ["proteomics-mcp@0.1.9", "serve"]
```

## Cursor

Add the server to `~/.cursor/mcp.json` (all projects) or `.cursor/mcp.json` (one project):

```json
{
  "mcpServers": {
    "proteomics-mcp": {
      "command": "uvx",
      "args": ["proteomics-mcp@0.1.9", "serve"]
    }
  }
}
```

## Windsurf

Add the same `mcpServers` entry to `~/.codeium/windsurf/mcp_config.json`, then refresh the MCP servers in Windsurf's settings.

## Other clients

Any MCP client that can start a stdio server works with the command above. [`install.json`](https://github.com/UadKLab/proteomics-mcp/blob/main/install.json) holds the same `mcpServers` entry for clients that import one.

## Settings

Set these as environment variables in the client's server entry (`"env": {...}` in JSON, an `env` table in TOML):

| Variable | Effect |
|---|---|
| `PROTEOMICS_MCP_OUTPUT_DIR` | Directory for relative output paths |
| `PROTEOMICS_MCP_LOG_LEVEL` | `DEBUG`, `INFO` (default), `WARNING` or `ERROR`; logs go to the client's MCP log |

## Check the connection

Ask the agent to list the proteomics-mcp tools; it should find 51. If the server does not start, the client's MCP log shows the error; see also the [troubleshooting prompts](prompt_examples/08_troubleshooting.md).
