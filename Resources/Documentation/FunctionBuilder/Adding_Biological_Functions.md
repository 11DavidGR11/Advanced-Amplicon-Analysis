# Adding biological functions

The Biological Function Builder creates declarative definitions in `Resources/FunctionalDB/CustomFunctions/`. Custom functions are displayed together in the **Custom functions** interface category, while their original scientific category is preserved as `biological_category`.

## Worked example

Create a function for nitrite respiration with either nitrite reductase plus nitric-oxide reductase.

| Field | Example value |
|---|---|
| Stable function ID | `nitrite_respiration_custom` |
| Display name | `Nitrite respiration` |
| Biological category | `Nitrogen metabolism` |
| Scientific description | `Potential reduction of nitrite to nitric oxide.` |
| Diagnostic genes | `nirK, nirS` |
| Supporting genes | `norB` |
| Accessory genes | `nosZ` |
| Decision rule | `Complex Boolean expression` |
| Boolean rule | `(nirK OR nirS) AND norB` |
| Interpretation | `Reports genomic potential, not activity.` |

The genes already exist in the standard Gene Dictionary, so the synonym field can remain empty. The saved function will appear under **Custom functions**, with **Nitrogen metabolism** retained as its biological category.

## Adding genes and synonyms

The Builder validates every declared marker against the combined built-in and custom Gene Dictionary. For a new canonical gene, add at least one accepted synonym using one line per gene:

```text
newGene = accepted product name; alternative annotation; historical name
otherGene = curated annotation name
```

Existing genes may also receive additional synonyms. Saving creates or updates:

```text
Resources/FunctionalDB/CustomGeneAliases.json
```

The previous file is backed up as `.bak`. A synonym cannot belong to two different canonical genes. Avoid broad terms that could match unrelated products.

## Rules

The Builder supports:

- **All diagnostic genes**: every diagnostic marker is required.
- **Minimum diagnostic genes**: at least the selected number of diagnostic markers is required.
- **Minimum genes across all roles**: counts diagnostic, supporting and accessory markers together.
- **Complex Boolean expression**: combines declared genes with `AND`, `OR`, `NOT` and parentheses.

Examples:

```text
(nirK OR nirS) AND norB
(geneA AND geneB) OR (geneC AND NOT geneD)
geneA AND (geneB OR geneC)
```

Only genes declared in the diagnostic, supporting or accessory fields may appear in the expression. Operator precedence is `NOT`, then `AND`, then `OR`; parentheses are recommended whenever the biological logic could be ambiguous.

## Validation and activation

Before saving, the Builder checks IDs, marker roles, thresholds, Boolean syntax, dictionary membership and synonym conflicts. Built-in function IDs cannot be overwritten. Existing custom definitions and dictionaries are backed up before replacement. Restart Biological Analysis after saving so the registry and dictionary are reloaded.

A positive result indicates putative genomic potential, not expression or activity. Complex rules still require scientific review, reference cases and appropriate QA.
