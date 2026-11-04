#!/bin/bash
# H7, Notification: a macOS notification when a session waits for the user.
msg=$(/usr/bin/python3 -c 'import json,sys; print(json.load(sys.stdin).get("message","Claude Code needs you")[:200])' 2>/dev/null)
/usr/bin/osascript -e "display notification \"${msg//\"/\'}\" with title \"Claude Code\"" >/dev/null 2>&1 || true
