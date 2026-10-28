# Usage with agentic workflows

minelastic can be run by an AI coding agent such as Claude Code. For this, the repository contains an agent
skill in the Agent Skills format,
[`skills/minelastic`](https://github.com/utksi/minelastic/tree/main/skills/minelastic): a `SKILL.md` file
with the rules and the commands, and two reference files that the agent opens when it needs them.

With the skill, an agent

* uses `minelastic.calculate` with ASE calculators, and `info`, `generate` and `collect` with DFT codes;
* chooses the options for the quantity you ask for, for example `--ionic` and `LCALCPOL` for the
  piezoelectric tensors with VASP;
* relaxes the reference structure first, uses the same settings for every cell, points out a cutoff or
  k-point mesh that looks too loose, and checks the symmetry before it writes any input;
* writes the inputs of a DFT code but runs the code only if you ask it to.

## Installation

In Claude Code:

<!-- skip: Claude Code commands -->
```text
/plugin marketplace add utksi/minelastic
/plugin install minelastic@minelastic
```

For other agents, copy the folder `skills/minelastic` into the skills folder of the agent, or point the agent
to `skills/minelastic/SKILL.md`.

## Example requests

<!-- input -->
```text
Compute the relaxed-ion elastic constants of POSCAR with MACE.

vasp_aln/ is my VASP template for wurtzite AlN. Write the inputs for the elastic constants and the
piezoelectric tensors into el/.

The runs in el/ have finished. Collect the elastic constants and check them.
```

The agent reports the constants in GPa together with the space group, the number of calculations, the checks
that minelastic prints and any assumption it made.

!!! note

    The agent needs minelastic installed ([Installation](../get-started/installation.md)). Check the inputs
    it writes as you would check your own; the skill does not replace a convergence study.
