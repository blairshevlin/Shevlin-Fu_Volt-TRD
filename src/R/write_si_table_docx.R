# Copyright (C) 2026 Blair Shevlin <blairshevlin@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 2026/11/16    Blair Shevlin                          optional .docx export of model tables, in the table format of the
#                                                      former generate_tablesSupplement_ntRegressions.R
#
# Used by main_analyses.R and supplementary_information.R only when export_docx = TRUE.
# Writes one Word file per table to results/si_tables/, ready to copy into the Supplementary Information.

write_si_table_docx = function(table, title, file_name) {
  out_dir = here::here("results", "si_tables")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  n_col = ncol(table)
  ft = flextable::flextable(as.data.frame(table)) |>
    flextable::theme_vanilla() |>
    flextable::bold(part = "header") |>
    flextable::bg(bg = "#f0f0f0", part = "header") |>
    flextable::align(align = "center", part = "header") |>
    flextable::align(j = 3:n_col, align = "center", part = "body") |>
    flextable::align(j = 1:2, align = "left", part = "body") |>
    flextable::fontsize(size = 9) |>
    flextable::font(fontname = "Times New Roman") |>
    flextable::width(j = 1, width = 1.5) |>
    flextable::width(j = 2, width = 2) |>
    flextable::width(j = 3:n_col, width = 0.7) |>
    flextable::valign(valign = "top") |>
    flextable::padding(padding.top = 2, padding.bottom = 2)

  officer::read_docx() |>
    officer::body_add_par(title, style = "heading 1") |>
    flextable::body_add_flextable(ft) |>
    print(target = file.path(out_dir, file_name))
  cat("Wrote", file.path(out_dir, file_name), "\n")
}
