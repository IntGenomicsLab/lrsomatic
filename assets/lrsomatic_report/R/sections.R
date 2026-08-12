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
