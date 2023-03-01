# Snekmer KBase Narrative App

# [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7671938.svg)](https://doi.org/10.5281/zenodo.7671938)

This [KBase](https://kbase.us) module was generated leveraging the KBase Software Development Kit (SDK) for running Snekmer modeling and search sequence anlaysis applications from in the KBase Narrative interface. You can also learn more about the apps implemented in this module from the [catalog page](https://narrative.kbase.us/#catalog/modules/Snekmer) or [spec file]($module_name.spec).

# Dependencies

Requires installation of the [KBase Software Development Kit (SDK)](https://github.com/kbase/kb_sdk) to use this module. May requires knowledgeable usage and workflow literacy of Snekmer software package listed under the citation reference listed below.

# Setup and test

Add your KBase developer token to `test_local/test.cfg` and run the following:

```bash
$ make
$ kb-sdk test
```

After making any additional changes to this repo, run `kb-sdk test` again to verify that everything still works.

# Installation from another module

To use this code in another SDK module, call `kb-sdk install Snekmer` in the other module's root directory.

# Help

To learn more about the KBase App Catalog visit the [KBase SDK User Documentation](https://kbase.github.io/kb_sdk_docs/) landing page. Additional [FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html) and supporting documentation can be found in the [KBase Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html).

# Citation Guidance

1. Jerger, Abby, Jason, McDermott E., & Nelson, William B. (2023). PerCon SFA Snekmer KBase Narrative App (v1.0.0). Zenodo. [https://doi.org/10.5281/zenodo.7671938](https://doi.org/10.5281/zenodo.7671938)
2. McDermott, Jason E., Chang, Christine H., Jerger, Abby, Nelson, William B., & Jacobson, Jeremy R. (2023). Snekmer: A scalable pipeline for protein sequence fingerprinting using amino acid recoding (AAR) (v1.0.3). Zenodo. [https://doi.org/10.5281/zenodo.7662597](https://doi.org/10.5281/zenodo.7662597)
