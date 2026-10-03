# Analysis code

This repository contains reviewed Python and R analysis scripts. Research data, generated results, manuscripts, and local configuration are excluded from version control.

Python preprocessing scripts optionally read `private_inputs.json` in the project root for `nmr_file`, `lipid_file`, `ethanol_file`, `dedicated_sample_prefix`, and `dedicated_basename`. Input filenames are relative to `raw_data/`. This configuration must remain local. Without it, generic filenames and sample prefixes are used.

Scripts require locally supplied inputs in `raw_data/`, `clean_data/`, and `outputs/`. Python preprocessing uses openpyxl and xlrd; R dependencies are listed in each script. The dementia analysis also requires a local report template.

Keep the explicit allowlist in `.gitignore`. Review source for credentials, embedded participant data, and confidential findings before adding new files.
