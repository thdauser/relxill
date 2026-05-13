# relxill — Developer Documentation

This site collects narrative documentation for the
[**relxill**](https://github.com/thdauser/relxill) model: notes on conventions,
physics derivations, normalization factors, and design decisions that go
beyond what fits in a single function header.

For build instructions, public usage, and citation policy see the
top-level `README.md` and `NAME_AND_CITATION_POLICY.md` in the repository.

## Topic Pages

- @subpage returnrad_normalizations — flux normalizations in the
  return-radiation reflection path
  (`calc_zone_direct_plus_xillver_reflected_flux`).

## API Reference

The full API reference is available via the sidebar
("Files", "Classes", "Namespaces"). Only entities with a Doxygen comment
are documented — `EXTRACT_ALL` is intentionally `NO` so this site stays
narrative-first and is not flooded with auto-generated stubs.

## Adding a Topic Page

1. Create `doc/<topic>.md` starting with `# Title {#topic_anchor}`.
2. Add a `- @subpage topic_anchor` line above.
3. Reference the page from the relevant function/file with
   `@see @ref topic_anchor`.

Math is rendered with MathJax: use `\f$ ... \f$` inline and
`\f[ ... \f]` for display equations.
