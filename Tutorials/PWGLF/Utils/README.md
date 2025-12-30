# PWGLF Utils Tutorial

This directory contains tutorials for common PWGLF utilities and helper classes.

## Collision Cuts Tutorial

**File:** `collisionCuts.cxx`

**Workflow:** `o2-analysis-collision-cuts-tutorial`

This tutorial demonstrates how to use the `collisionCutsGroup` pattern for automatic configurable registration of collision event selection parameters.

### Key Features

- **Automatic Configurable Registration**: All collision selection parameters are automatically registered by the DPL framework with a single line of code
- **Reduced Boilerplate**: Traditional approach requires ~30 lines of manual configurable declarations, the new pattern reduces this to ~5 lines
- **Comprehensive Event Selection**: Includes z-vertex cuts, centrality cuts, time frame cuts, occupancy cuts, and various Run 3 event selection flags

### What You'll Learn

1. How to declare a `collisionCutsGroup` configurable
2. How to apply settings from the group to a `CollisonCuts` object
3. How to use the collision cuts checker to filter events
4. How to fill QA histograms for selected events

### Customizing Cuts

The collision cuts can be customized via JSON configuration file:

```json
{
  "collision-cuts-tutorial": {
    "collCuts": {
      "VtxRange": 10.0,
      "CentMin": 0.0,
      "CentMax": 90.0,
      "CentEstimator": 1,
      "IsRun3": true,
      "DoOccupancySel": true,
      "OccupancyMax": 1000
    }
  }
}
```

For a complete list of available parameters, run:
```bash
o2-analysis-collision-cuts-tutorial --help full
```

### Related Documentation

- Main collision cuts implementation: `PWGLF/Utils/collisionCuts.h`
- Configurable group implementation: `PWGLF/Utils/collisionCutsGroup.h`
- Event selection parameters: `Common/CCDB/EventSelectionParams.h`
