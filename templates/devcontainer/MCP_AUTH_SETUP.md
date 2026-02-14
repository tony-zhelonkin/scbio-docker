# Claude Code + MCP: Complete Setup & Security Guide

## Table of Contents
1. [Quick Start](#quick-start)
2. [Background: How MCP Works](#background-how-mcp-works)
3. [Security Implications & Risk Analysis](#security-implications--risk-analysis)
4. [Attack Vectors & Mitigations](#attack-vectors--mitigations)
5. [Implementation Routes](#implementation-routes)
6. [Automated Setup (Recommended)](#automated-setup-recommended)
7. [Manual OAuth Authentication](#manual-oauth-authentication)
8. [Troubleshooting](#troubleshooting)
9. [References](#references)

---

## Quick Start

### Automated Setup (Already Configured!)

This repository includes automated Claude Code + MCP setup. When you rebuild/restart the container:

1. **Claude Code** is automatically installed
2. **MCP dependencies** (`uv`, Node.js) are automatically installed
3. **MCP servers** are configured from `.mcp.json`
4. All you need to do manually is **authenticate context7** (OAuth)

### Next Steps After Container Start

1. **Set up SSH port forwarding** (one-time, see [Manual OAuth Authentication](#manual-oauth-authentication)):
```bash
# On your LOCAL machine
ssh -L 45454:localhost:45454 your-username@your-remote-server
```

2. **Start Claude Code** inside the container:
```bash
claude
```

1. **Authenticate context7** (inside Claude Code):
```
/mcp
```
   - Follow prompts to authenticate context7 via browser
   - Other MCP servers (serena, sequential-thinking) work without authentication

---

## Background: How MCP Works

### What is Model Context Protocol (MCP)?

**MCP** is a protocol that allows AI assistants like Claude to interact with external tools and data sources through a standardized interface. Think of it as "plugins for AI agents."

```
┌─────────────┐         ┌─────────────┐         ┌──────────────────┐
│             │         │             │         │                  │
│  Claude AI  │ ◄─────► │ MCP Client  │ ◄─────► │  MCP Server(s)   │
│  (LLM)      │         │ (Claude     │         │  (Tools/Data)    │
│             │         │  Code)      │         │                  │
└─────────────┘         └─────────────┘         └──────────────────┘
                                                         │
                                                         ├─► context7
                                                         ├─► serena
                                                         └─► sequential-thinking
```

### MCP Transport Mechanisms

MCP servers communicate via three transport types:

#### 1. **stdio** (Standard Input/Output)
- **How it works**: MCP server runs as a subprocess, communicates via stdin/stdout
- **Security**: Runs with container user's permissions, local-only access
- **Examples**: `serena`, `sequential-thinking`
- **Pros**: Simple, fast, no network exposure
- **Cons**: Can't communicate across containers/hosts without proxying

```javascript
// Example stdio configuration
{
  "type": "stdio",
  "command": "uvx",
  "args": ["--from", "git+https://github.com/oraios/serena", "serena", "start-mcp-server"]
}
```

#### 2. **SSE** (Server-Sent Events) - **DEPRECATED**
- **How it works**: HTTP endpoint streaming events to client
- **Security**: Network-exposed port, requires authentication
- **Examples**: `context7` (older implementation)
- **Pros**: Can work across network boundaries
- **Cons**: Being replaced by Streamable HTTP, security risks if unauthenticated

#### 3. **Streamable HTTP** (Modern, replacing SSE)
- **How it works**: Bi-directional HTTP streaming
- **Security**: Same as SSE but more robust protocol
- **Examples**: Newer MCP servers
- **Pros**: Better error handling, standardized protocol
- **Cons**: Still requires network security considerations

### How MCP Servers Access Your Code

When you enable an MCP server, Claude Code grants it specific **tools** that can:

1. **Read files**: Access project files, source code, documentation
2. **Search code**: Grep, semantic search, AST parsing
3. **Execute commands**: Run shell commands (depending on server)
4. **Network requests**: Fetch external documentation, APIs
5. **Modify files**: Some MCP servers can edit/create files

**Critical Understanding**: MCP servers have the **same permissions as the user running Claude Code**. In a container, this means:
- Can read/write all files the container user can access
- Can execute any command the container user can run
- Can make network requests from the container
- **Cannot** escape the container (unless there's a container vulnerability)

---

## Security Implications & Risk Analysis

### Core Security Principle

> **MCP servers are third-party code running with your permissions.**

This is fundamentally similar to:
- Installing a browser extension
- Running `npm install` from an unknown package
- Adding a Python library via `pip install`
- Executing a bash script from the internet

### What Can MCP Servers Do?

#### 1. **Code & Data Access**
MCP servers can access:
- ✅ **Your entire codebase** (all files in `/workspaces/DC_Dictionary`)
- ✅ **Environment variables** (API keys, secrets if exposed)
- ✅ **Git history** (including deleted secrets)
- ✅ **Mounted volumes** (read-only data directories in your setup)
- ❌ **Files outside container** (protected by Docker isolation)

**Risk Level**: 🟡 **MEDIUM** - Limited to container, but full project access

**Mitigation**:
- Never commit secrets to git (use `.gitignore`, environment files)
- Review mounted volumes in `docker-compose.yml`
- Use read-only mounts for sensitive data (already done: `:ro` flag)

#### 2. **Network Access**
MCP servers can:
- ✅ **Fetch external URLs** (documentation, APIs)
- ✅ **Exfiltrate data** (send code/secrets to external servers)
- ✅ **Download malicious payloads**
- ⚠️ **Access internal network** (if container networking allows)

**Risk Level**: 🔴 **HIGH** - Data exfiltration is the primary risk

**Mitigation**:
- Only use **trusted MCP servers** (official/well-reviewed)
- Monitor network traffic (use `docker stats`, firewall rules)
- Use network policies to restrict container egress (advanced)
- Review MCP server source code before use

#### 3. **Execution Capabilities**
Some MCP servers can:
- ✅ **Run shell commands** (arbitrary code execution)
- ✅ **Modify files** (delete, edit, create)
- ✅ **Install packages** (npm, pip, etc.)
- ❌ **Gain root** (runs as devuser, not root)

**Risk Level**: 🔴 **HIGH** - Code execution within container

**Mitigation**:
- Use stdio-based MCP servers (local execution only)
- Review MCP server capabilities before enabling
- Keep container user non-root (already done: `user: "${LOCAL_UID}:..."`)
- Use separate containers for untrusted MCP servers

#### 4. **OAuth & Authentication Tokens**
OAuth-based MCP servers (like context7):
- ✅ **Store tokens** in `~/.claude.json` (container-local)
- ⚠️ **Tokens valid across sessions** (until revoked)
- ⚠️ **Token theft via container compromise** (if attacker gains shell access)

**Risk Level**: 🟡 **MEDIUM** - Tokens are scoped to MCP services

**Mitigation**:
- Tokens stored in `~/.claude.json` (not version-controlled)
- Revoke tokens via `/mcp` → "Clear authentication"
- Use SSH port forwarding (not port exposure) for OAuth callback
- Regularly audit authorized applications in OAuth provider

---

## Attack Vectors & Mitigations

### Attack Vector 1: Malicious MCP Server

**Scenario**: You install an MCP server from an untrusted source that:
1. Reads your entire codebase
2. Scans for API keys, credentials, secrets
3. Exfiltrates data to attacker's server
4. Injects backdoor code into your project

**Likelihood**: 🟡 **MEDIUM** (depends on source trustworthiness)

**Mitigations**:
- ✅ **Only use official/vetted MCP servers**:
  - `context7`: Official, hosted by context7.com
  - `serena`: Open-source, GitHub repo reviewable
  - `sequential-thinking`: Official Anthropic/MCP server
- ✅ **Review source code** before installation (especially stdio servers)
- ✅ **Check repository stars/activity** (community trust signals)
- ✅ **Use network monitoring** to detect unexpected outbound connections
- ✅ **Run untrusted MCP servers in separate containers** (advanced)

**Detection**:
```bash
# Monitor container network activity
docker stats dev-core --no-stream
# Check active network connections
sudo netstat -tulpn | grep $(docker inspect -f '{{.State.Pid}}' <container-id>)
```

### Attack Vector 2: Supply Chain Attack

**Scenario**: A legitimate MCP server's dependencies are compromised:
1. MCP server uses `npm`/`pip` package with malicious code
2. Package update introduces backdoor
3. Your container installs compromised version
4. Backdoor activates on MCP server start

**Likelihood**: 🟡 **MEDIUM** (same as any npm/pip package)

**Mitigations**:
- ✅ **Pin versions** in `.mcp.json` (when possible)
- ✅ **Review dependency trees** (`npm ls`, `pip show`)
- ✅ **Use lockfiles** for reproducible builds
- ✅ **Scan for vulnerabilities** (`npm audit`, `pip-audit`)
- ✅ **Container rebuilds get fresh dependencies** (this setup does this)

**Current Setup**: Our automated script re-installs MCP dependencies on every container start, which:
- ✅ **Pros**: Always get latest security patches
- ❌ **Cons**: Could pull compromised version if attack happens

**Considering**: to pin versions in `setup_claude_mcp.sh`

### Attack Vector 3: OAuth Token Theft

**Scenario**: Attacker gains access to container and steals OAuth tokens:
1. Attacker exploits vulnerability to gain shell access
2. Reads `~/.claude.json` (contains OAuth tokens)
3. Uses stolen tokens to access context7 with your identity
4. Queries sensitive documentation, exfiltrates data

**Likelihood**: 🟢 **LOW** (requires prior container compromise)

**Mitigations**:
- ✅ **Token file permissions**: `~/.claude.json` is `0600` (user-only read)
- ✅ **Non-persistent home directory**: Tokens lost on container rebuild
- ✅ **Regular token revocation**: `/mcp` → "Clear authentication"
- ✅ **SSH port forwarding** (not exposed port): Callback only via SSH tunnel

**Unique to This Setup**: Container rebuilds automatically **invalidate all tokens** (home directory is not volume-mounted), forcing re-authentication. This is a **security feature**.

### Attack Vector 4: Man-in-the-Middle (OAuth Callback)

**Scenario**: Attacker intercepts OAuth callback to steal authorization code:
1. You authenticate context7 via browser
2. OAuth callback redirects to `http://localhost:45454/callback?code=...`
3. Attacker intercepts traffic (if callback goes over public network)
4. Steals authorization code, exchanges for access token

**Likelihood**: 🟢 **LOW** (mitigated by SSH tunnel)

**Mitigations**:
- ✅ **SSH port forwarding**: Callback traffic encrypted via SSH tunnel
- ✅ **Local browser**: OAuth happens on local machine, not remote
- ✅ **Short-lived authorization codes**: Codes expire in seconds
- ❌ **Don't expose port 45454 publicly**: Never do `0.0.0.0:45454` on remote server

**Current Setup**: OAuth callback flow:
```
Local Browser
    │
    └─► https://mcp.context7.com/auth (HTTPS, encrypted)
         │
         └─► callback: http://localhost:45454/... (local-only)
              │
              └─► SSH tunnel (encrypted) → Remote container
```

All traffic is encrypted except the local `localhost` portion (which never leaves your machine).

### Attack Vector 5: Compromised Workstation

**Scenario**: Your local machine (not the container) is compromised:
1. Malware on your laptop steals SSH keys
2. Attacker gains access to remote server via stolen keys
3. Attacker accesses dev container, steals code/data
4. MCP servers are irrelevant (attack bypasses them)

**Likelihood**: 🟡 **MEDIUM** (depends on local security practices)

**Mitigations**:
- ✅ **SSH key passphrase**: Always encrypt SSH keys
- ✅ **2FA for SSH**: Use hardware keys (YubiKey) or OTP
- ✅ **Host-based firewall**: Limit SSH access by IP
- ✅ **Regular security updates**: Keep local OS patched
- ✅ **Separate work/personal machines**: Don't mix contexts

**MCP-Specific**: If your workstation is compromised, MCP security is moot. Focus on:
1. Securing your local machine first
2. Using SSH best practices (keys, 2FA, known_hosts)
3. Monitoring SSH access logs on remote server

### Attack Vector 6: Container Escape

**Scenario**: MCP server exploits Docker vulnerability to escape container:
1. MCP server triggers container runtime bug
2. Gains access to host system (remote server)
3. Accesses other containers, host data
4. Installs persistent backdoor on host

**Likelihood**: 🟢 **VERY LOW** (Docker is well-hardened)

**Mitigations**:
- ✅ **Run as non-root**: Container user is `devuser` (UID 1000)
- ✅ **Read-only mounts**: Data volumes mounted `:ro`
- ✅ **Keep Docker updated**: Use latest Docker version
- ✅ **Resource limits**: `cpus`, `memory` limits prevent DoS
- ⚠️ **No AppArmor/SELinux profile**: Considering to adding for paranoia

**Current Setup**: Standard Docker security, no custom hardening

---

## Implementation Routes

### Route A: postStartCommand (RECOMMENDED) ✅ **IMPLEMENTED**

**What it does**:
- Runs `.devcontainer/scripts/setup_claude_mcp.sh` every time container starts
- Installs Claude Code, uv, Node.js if not present
- Creates `.mcp.json` if not present
- Idempotent (safe to run multiple times)

**Pros**:
- ✅ Fully automated, zero manual steps (except OAuth)
- ✅ Always uses latest versions (auto-updates)
- ✅ Works on any machine (no pre-built images needed)
- ✅ Non-persistent home directory (security: tokens don't survive rebuilds)
- ✅ Easy to audit (single script)

**Cons**:
- ❌ Slight startup delay (~30s first time, ~5s subsequent)
- ❌ Requires internet (downloads dependencies)
- ❌ No version pinning (could break if upstream changes)

**When to use**: Default for most users, development environments

**This is what we've implemented!**

### Route B: postStartCommand + Volume Mount

**What it does**:
- Same as Route A, but adds volume mount for `~/.local` and `~/.claude`
- Persists Claude Code installation and OAuth tokens across rebuilds

**Implementation**:
```yaml
# In docker-compose.yml
volumes:
  - claude-home:/home/devuser/.local
  - claude-config:/home/devuser/.claude

volumes:
  claude-home:
  claude-config:
```

**Pros**:
- ✅ Faster startup (no re-installation)
- ✅ OAuth tokens persist (don't need to re-authenticate)

**Cons**:
- ❌ **Security risk**: Tokens survive container compromise
- ❌ **Stale versions**: Manual updates required
- ❌ **Shared state**: Tokens/config shared across container rebuilds

**When to use**: Production-like environments where startup speed matters more than security

**Security trade-off**: Persistence = convenience but **reduced security**

### Route C: Modify Base Image (Advanced)

**What it does**:
- Rebuild `scdock-r-dev:v0.5.1` image with Claude Code baked in
- Create `Dockerfile` extending base image

**Implementation**:
```dockerfile
# .devcontainer/Dockerfile
FROM scdock-r-dev:v0.5.1

# Install Claude Code dependencies
RUN apt-get update && apt-get install -y curl gnupg && \
    curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    pip3 install uv

# Install Claude Code
RUN curl -fsSL https://claude.ai/install.sh | bash -s latest

USER devuser
```

**Pros**:
- ✅ Fastest startup (everything pre-installed)
- ✅ Version controlled (Dockerfile in git)
- ✅ Reproducible builds

**Cons**:
- ❌ Requires image rebuild on updates
- ❌ More complex CI/CD pipeline
- ❌ Larger image size
- ❌ Not portable (requires custom image registry)

**When to use**: Team environments with shared infrastructure

---

## Automated Setup (Recommended)

### How It Works

This repository uses **Route A** (postStartCommand). Here's what happens:

1. **Container starts** (rebuild or restart)
2. **`postStartCommand` triggers** in `devcontainer.json`:
   ```json
   "postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
   ```
3. **Setup script runs** (`.devcontainer/scripts/setup_claude_mcp.sh`):
   - Checks if Claude Code installed → installs if missing
   - Checks if `uvx` installed → installs `uv` if missing
   - Checks if `npx` installed → installs Node.js if missing
   - Creates `.mcp.json` if missing
   - Runs sanity checks
4. **Sanity checks run** (`.devcontainer/scripts/poststart_sanity.sh`)
5. **Container ready** with Claude Code + MCP servers configured

### What You Need to Do

**One-time setup** (first container start):

1. **Wait for setup script to complete** (watch terminal output)
2. **Set up SSH port forwarding**:
   ```bash
   # On your LOCAL machine
   ssh -L 45454:localhost:45454 your-username@your-remote-server
   ```
3. **Start Claude Code**:
   ```bash
   claude
   ```
4. **Authenticate context7**:
   ```
   /mcp
   ```
   - Select "context7"
   - Click the OAuth URL in your local browser
   - Authorize the application
   - Context7 authenticated ✅

**Subsequent container starts**:
- SSH port forwarding still needed (or add to `~/.ssh/config`)
- OAuth tokens expire → re-authenticate when prompted
- Everything else is automatic

---

## Manual OAuth Authentication

### Why Manual Steps Are Needed

OAuth authentication for `context7` requires:
1. Opening a browser (not available in SSH/container)
2. Callback URL accessible from browser (`http://localhost:45454/callback`)
3. Container port accessible from your local machine

**Without port forwarding**: Browser can't reach container's port 45454
**With port forwarding**: Local browser → SSH tunnel → Container ✅

---

## Port Forwarding 101: Complete Beginner's Guide

### What is Port Forwarding?

Think of port forwarding like creating a **secret tunnel** between two computers. In our case:
- **Your laptop** (where your browser runs) needs to talk to
- **The dev container** (running on a remote server via SSH)

Normally, your browser can only access websites on the internet OR things running on your laptop (`localhost`). It **cannot** reach into remote servers or containers. Port forwarding solves this by creating a tunnel.

### Visual Explanation

**WITHOUT Port Forwarding** (doesn't work):
```
┌─────────────────┐                    ┌──────────────────┐
│  Your Laptop    │                    │  Remote Server   │
│                 │                    │                  │
│  🌐 Browser     │   ❌ Can't reach   │  🐳 Container    │
│                 │   ════════════════▶│  (port 45454)   │
│                 │                    │                  │
└─────────────────┘                    └──────────────────┘
    (local)                                 (remote)
```

**WITH Port Forwarding** (works!):
```
┌─────────────────┐         SSH Tunnel (encrypted)        ┌──────────────────┐
│  Your Laptop    │         ┌───────────────────┐         │  Remote Server   │
│                 │         │                   │         │                  │
│  🌐 Browser ────┼────────▶│  localhost:45454  │────────▶│  🐳 Container    │
│                 │         │  (your laptop)    │         │  (port 45454)   │
│                 │         └───────────────────┘         │                  │
└─────────────────┘                                       └──────────────────┘
    (local)                                                    (remote)

    Your browser thinks it's talking to localhost,
    but SSH secretly forwards the traffic to the remote container!
```

### The Magic: What `-L 45454:localhost:45454` Means

The SSH command `ssh -L 45454:localhost:45454 user@remote` has three parts:

```
ssh -L 45454:localhost:45454 user@remote
       ↑     ↑       ↑
       │     │       │
       │     │       └─ Port on remote server/container
       │     └───────── Where traffic goes (remote's localhost)
       └─────────────── Port on YOUR laptop
```

**In plain English**:
> "SSH, create a tunnel so that when something on MY laptop tries to connect to `localhost:45454`, secretly send that traffic through the SSH connection to the remote server's `localhost:45454`, which Docker then forwards to the container's port 45454."

### Real-World Flow for OAuth Authentication

Here's **exactly** what happens when you authenticate context7:

1. **You run** (on remote server, inside container):
```bash
claude
/mcp
# Select "context7" → OAuth URL appears
```

2. **Claude Code says**: "Go to this URL: `https://mcp.context7.com/auth?callback=http://localhost:45454/callback`"

3. **You copy URL and paste in YOUR LAPTOP's browser**

4. **Browser opens** `https://mcp.context7.com/auth?...`
   - You log in to context7
   - You authorize the application
   - context7 says: "OK, redirecting you to `http://localhost:45454/callback?code=SECRET123`"

5. **Browser tries to connect to** `http://localhost:45454/callback?code=SECRET123`
   - Normally this would fail (nothing running on your laptop's port 45454)
   - **BUT** SSH port forwarding intercepts this!
   - SSH says: "Aha! I know where port 45454 really is!"

6. **SSH tunnel forwards the request**:
   - Traffic goes through encrypted SSH connection
   - Arrives at remote server's `localhost:45454`
   - Docker forwards to container's port 45454
   - Claude Code (listening on port 45454 inside container) receives it!

7. **Claude Code gets the authorization code** and completes authentication ✅

### Step-by-Step: SSH Port Forwarding

#### Prerequisites Check

Before starting, verify you know:
1. **Your remote server's address**: e.g., `10.1.2.3` or `my-server.university.edu`
2. **Your username on remote server**: e.g., `user` or `devuser`
3. **You can SSH to remote server**: `ssh username@remote-server` works

#### Option 1: Command-Line Port Forwarding (Per Session Quick Test)

On your **local machine**, connect to remote server with port forwarding:

```bash
ssh -L 45454:localhost:45454 your-username@your-remote-server
```

**What this does**:
- Forwards local port 45454 → remote server's localhost:45454
- Remote server's localhost:45454 → container port 45454 (via docker-compose mapping)
- Creates encrypted tunnel for OAuth callback

**Keep this SSH session open** while using Claude Code.

#### Option 2: SSH Config (Persistent)

On your **local machine**, edit `~/.ssh/config`:

```bash
Host your-remote-server
    HostName your.server.address.com
    User your-username
    LocalForward 45454 localhost:45454
```

Now **every SSH connection** automatically forwards port 45454:

```bash
ssh your-remote-server  # Port forwarding happens automatically
```

---

### Complete Worked Example (Step-by-Step)

Let's walk through a **real example** from start to finish. This is what you'll actually type and see:

#### Example Setup:
- **Your laptop**: MacBook (or Windows/Linux)
- **Remote server**: `research-server.university.edu`
- **Your username**: `user`
- **Container name**: `dev-core` (from docker-compose.yml)

#### Step 1: Open Terminal on Your Laptop

```bash
# On YOUR LAPTOP (Terminal or PowerShell)
# Check you can reach the remote server first
ssh user@research-server.university.edu

# You should get a shell on the remote server
# Type 'exit' to close and come back to your laptop
exit
```

#### Step 2: Connect with Port Forwarding

```bash
# On YOUR LAPTOP
# This time, add the -L flag for port forwarding
ssh -L 45454:localhost:45454 user@research-server.university.edu

# You'll see your normal SSH login
# No visual difference, but port 45454 is now tunneled!
```

**What you'll see**:
```
Welcome to research-server.university.edu
Last login: Mon Oct 21 10:00:00 2025 from 10.1.2.3

user@research-server:~$
```

**⚠️ IMPORTANT**: **Keep this terminal window open!** The port forwarding only works while the SSH session is active.

#### Step 3: Attach to the Container

```bash
# On REMOTE SERVER (in the SSH session from step 2)
# Navigate to your project
cd /path/to/DC_Dictionary

# Attach to the running container
docker exec -it dev-core bash

# You're now INSIDE the container
devuser@95319867f4a7:/workspaces/DC_Dictionary$
```

#### Step 4: Start Claude Code

```bash
# INSIDE CONTAINER
# Start Claude Code
claude

# You'll see the Claude Code interface load
# Wait for it to be ready (shows "Claude Code >")
```

#### Step 5: Trigger OAuth Authentication

```
# Inside Claude Code prompt
/mcp
```

**What you'll see**:
```
Available MCP Servers:
  ⚪ context7 - Not authenticated
  ✅ serena - Ready
  ✅ sequential-thinking - Ready

Use arrow keys to select. Press Enter to authenticate.
```

Use arrow keys to select **context7**, press **Enter**.

#### Step 6: Copy OAuth URL

**Claude Code will display**:
```
To authenticate context7, open this URL in your browser:

https://mcp.context7.com/auth?client_id=...&redirect_uri=http://localhost:45454/callback

Waiting for authentication...
```

**Copy the entire URL** (Ctrl+C or Cmd+C).

#### Step 7: Open URL in YOUR LAPTOP's Browser

**⚠️ CRITICAL**: Open the browser **ON YOUR LAPTOP**, not on the remote server!

1. **Open browser** (Chrome, Firefox, Safari, whatever you use)
2. **Paste the URL** into the address bar
3. **Press Enter**

**What happens in browser**:
```
1. Browser opens: https://mcp.context7.com/auth?...
2. You see context7 login page
3. Log in with your credentials
4. You see: "Authorize Claude Code to access context7?"
5. Click "Allow" or "Authorize"
6. Browser redirects to: http://localhost:45454/callback?code=abc123...
```

**What you'll see**:
- Browser might show a loading page briefly
- Then displays: "Authentication successful! You can close this window."

#### Step 8: Verify Authentication

**Back in Claude Code** (in the SSH session):

You should see:
```
✅ Authentication successful!

context7 is now authenticated.
```

Type `/mcp` again to verify:

```
/mcp

Available MCP Servers:
  ✅ context7 - Authenticated ✅
  ✅ serena - Ready
  ✅ sequential-thinking - Ready
```

🎉 **Done!** All three MCP servers are now working!

#### What Just Happened?

Let's trace the **exact path** of the OAuth callback:

```
1. Browser (laptop) tries to connect to localhost:45454
                    ↓
2. SSH intercepts (because of -L 45454:localhost:45454)
                    ↓
3. SSH forwards through encrypted tunnel to remote server
                    ↓
4. Remote server's localhost:45454
                    ↓
5. Docker port mapping (docker-compose.yml ports: "45454:45454")
                    ↓
6. Container's port 45454
                    ↓
7. Claude Code (listening on port 45454 inside container)
                    ↓
8. ✅ Authentication complete!
```

---

### Common Mistakes & How to Avoid Them

#### ❌ Mistake 1: Forgetting Port Forwarding

**Wrong**:
```bash
# Just regular SSH (no -L flag)
ssh user@research-server.university.edu
```

**What happens**: OAuth callback fails with "Connection refused"

**Fix**: Use `-L 45454:localhost:45454`

#### ❌ Mistake 2: Opening OAuth URL on Remote Server

**Wrong**: Trying to open URL in a remote browser or using `curl`

**What happens**: No browser available, or callback doesn't reach Claude Code

**Fix**: **Always** open OAuth URL in your **laptop's browser**

#### ❌ Mistake 3: Closing SSH Session Too Early

**Wrong**:
```bash
ssh -L 45454:localhost:45454 user@remote
# Do some stuff
exit  # ← Closes SSH session BEFORE authenticating
```

**What happens**: Port forwarding dies, OAuth callback fails

**Fix**: Keep SSH session open until **after** authentication completes

#### ❌ Mistake 4: Using Wrong Port

**Wrong**:
```bash
ssh -L 8080:localhost:45454 user@remote  # Wrong local port
```

**What happens**: Browser tries localhost:45454, but SSH is listening on 8080

**Fix**: Use **same port number** on both sides: `45454:localhost:45454`

---

### How to Verify Port Forwarding is Working

**Before** trying OAuth, test if port forwarding is active:

#### On Your Laptop (New Terminal):

```bash
# While SSH session is still open, open a NEW terminal on your laptop
# Try to connect to localhost:45454
curl http://localhost:45454

# If port forwarding is working:
# - You should get SOME response (could be error from Claude Code, but connection succeeds)
# - If you get "Connection refused", port forwarding is NOT active
```

#### Inside the Container:

```bash
# Check if Claude Code is listening on port 45454
netstat -tuln | grep 45454

# Should show something like:
# tcp  0  0  0.0.0.0:45454  0.0.0.0:*  LISTEN
```

If both checks pass, port forwarding is configured correctly! ✅

---

### Step-by-Step: Authenticating context7

1. **Ensure SSH port forwarding is active** (see above)

2. **Start Claude Code** in container:
```bash
claude
```

1. **Open MCP menu**:
```
/mcp
```

1. **Authenticate context7**:
   - Use arrow keys to select "context7"
   - Press Enter
   - Claude Code displays OAuth URL (starts with `https://mcp.context7.com/auth?...`)

2. **Copy URL and open in LOCAL browser**:
   - **Do NOT open in remote server** (no GUI browser)
   - Open on your laptop/desktop
   - URL goes to context7, redirects to `http://localhost:45454/callback?code=...`
   - Callback reaches container via SSH tunnel ✅

3. **Authorize in browser**:
   - Follow context7's OAuth flow
   - Grant permissions
   - Browser redirects to `localhost:45454` (your local machine)
   - SSH tunnel forwards to remote container
   - Claude Code receives authorization ✅

4. **Verify authentication**:
```
/mcp
```
   - context7 should show "Authenticated ✅"

### Troubleshooting OAuth

**Error: "Connection refused" on callback**
- ✅ Check SSH port forwarding is active: `ssh -L 45454:localhost:45454 ...`
- ✅ Check container port mapping: `docker ps` should show `0.0.0.0:45454->45454/tcp`
- ✅ Check firewall allows port 45454

**Error: "Timeout waiting for callback"**
- ✅ Did you open URL in **local browser** (not remote)?
- ✅ Is SSH tunnel still active? (check terminal)
- ✅ Try reauthenticating: `/mcp` → "Clear authentication" → retry

**Error: "Invalid callback URL"**
- ✅ Check callback URL is `http://localhost:45454/...` (not `http://0.0.0.0:...`)
- ✅ Ensure you didn't modify port mapping to `0.0.0.0:45454`

---

## Troubleshooting

### Setup Script Fails

**Error: "Claude Code installation failed"**

Check:
```bash
# Ensure curl works
curl -I https://claude.ai/install.sh

# Check disk space
df -h

# Run installer manually
curl -fsSL https://claude.ai/install.sh | bash -s latest
```

**Error: "uv installation failed"**

Check:
```bash
# Ensure pip3 works
pip3 --version

# Install manually
sudo pip3 install uv

# Verify
which uvx
```

**Error: "Node.js installation failed"**

Check:
```bash
# Ensure NodeSource repo is accessible
curl -fsSL https://deb.nodesource.com/setup_20.x

# Install manually
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs

# Verify
which npx
```

### MCP Servers Not Working

**`serena` fails to start**

Check:
```bash
# Test manually
uvx --from git+https://github.com/oraios/serena serena --help

# Check Python version (needs 3.8+)
python3 --version

# Check git is available
git --version
```

**`sequential-thinking` fails to start**

Check:
```bash
# Test manually
npx -y @modelcontextprotocol/server-sequential-thinking

# Check npm registry access
npm config get registry

# Clear npm cache
npm cache clean --force
```

**`context7` authentication fails**

See [Manual OAuth Authentication](#manual-oauth-authentication) above.

### Container Performance Issues

**Slow startup (setup script takes minutes)**

This is normal on **first run** (downloads dependencies). Subsequent runs are faster due to caching.

To speed up:
- Use Route B (volume mounts) for persistent installation
- Or pre-build image with Route C (Dockerfile)

**High memory usage**

Check:
```bash
# Container stats
docker stats dev-core

# MCP server processes
ps aux | grep -E "(uvx|npx|node)"
```

MCP servers (especially serena) can be memory-intensive. Your container has 300GB limit, so this shouldn't be an issue.

---

## References

### Official Documentation

- **Claude Code Docs**: https://docs.claude.com/en/docs/claude-code
- **MCP Protocol Spec**: https://spec.modelcontextprotocol.io/
- **MCP Troubleshooting**: https://docs.claude.com/en/docs/claude-code/troubleshooting
- **Context7 MCP Server**: https://context7.com/
- **Serena MCP Server**: https://github.com/oraios/serena
- **Sequential Thinking MCP**: https://www.npmjs.com/package/@modelcontextprotocol/server-sequential-thinking

### Security Resources

- **Docker Security Best Practices**: https://docs.docker.com/engine/security/
- **SSH Port Forwarding Guide**: https://www.ssh.com/academy/ssh/tunneling/example
- **OAuth 2.0 Security**: https://oauth.net/2/
- **Supply Chain Security**: https://slsa.dev/

### Related Blog Posts

- **Self-hosted MCP in Dev Containers**: https://blairwadman.com/articles/self-hosted-mcp-server-in-development-containers/
- **Building AI-powered development containers**: https://medium.com/@brett_4870/building-a-secure-ai-development-environment-containerized-claude-code-mcp-integration-e2129fe3af5a
- **MCP Security Considerations**: (community resources)

### Community

- **Claude Code GitHub Issues**: https://github.com/anthropics/claude-code/issues
- **MCP GitHub**: https://github.com/modelcontextprotocol

---

**Questions or issues?** Consult `.devcontainer/scripts/setup_claude_mcp.sh` for implementation details.