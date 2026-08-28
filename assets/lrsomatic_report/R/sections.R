# Section-module contract: each report section registers a descriptor with
# id, title, locate(sample_dir, sample_id), and parse(inputs, section_data).
# See CLAUDE.md "Section-module contract" for the recipe to add a new section.

SECTIONS = list()

register_section = function(descriptor) {
  SECTIONS[[descriptor$id]] <<- descriptor
}

# Standard "nothing to show" notice used by section presentation shims.
section_notice = function(msg) {
  tags$div(class = "alert alert-info", msg)
}

# Provenance and caveats that belong under a table, folded away by default.
#
# This text is the only thing that surfaces a silently broken join — "0 of 49,003
# variants carry a VAF" is how the AF/VAF tag bug was caught — so it is collapsed rather
# than deleted. `summary` names what is inside, so a reader can tell whether it is worth
# opening without opening it.
table_details = function(..., summary = "Details") {
  tags$details(
    class = "table-details",
    tags$summary(summary),
    tags$div(class = "table-footnote", ...)
  )
}
