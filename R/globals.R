# Suppresses R CMD check NOTEs for data.table's non-standard evaluation.
# These symbols (column names and data.table's special variables .N/.SD/.)
# are resolved at runtime inside `[.data.table` calls; they are not actual
# unbound globals.
utils::globalVariables(c(
  ".", ".N", ".SD", "age", "id", "sex", "superpod"
))
