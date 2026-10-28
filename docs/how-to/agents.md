# AI coding agents

The repository contains an agent skill for minelastic in the Agent Skills format, in
[`skills/minelastic/`](https://github.com/utksi/minelastic/tree/main/skills/minelastic): `SKILL.md` with the
rules and the commands, and two reference files that the agent reads when it needs them. The skill tells
the agent when to use `calculate` and when the command line, which options give which quantity, and what to
check before and after a calculation: a stress-free reference, one set of settings for all cells, converged
settings and the symmetry.

## Claude Code

<!-- skip: Claude Code commands -->
```text
/plugin marketplace add utksi/minelastic
/plugin install minelastic@minelastic
```

The skill is then used whenever a request involves elastic constants, piezoelectric tensors or strained
inputs.

## Other agents

Copy the folder `skills/minelastic` into the skills folder of the agent, or point the agent at
`skills/minelastic/SKILL.md`.

!!! note

    The agent needs minelastic installed ([Installation](../get-started/installation.md)). The skill tells
    it not to run DFT codes unless you ask: it writes the inputs, and you run them as usual.
