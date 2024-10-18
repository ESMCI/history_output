set(HISTORY_HASH_SOURCES
  hist_hash_table.F90
  hist_hashable.F90
)

target_sources(history PRIVATE ${HISTORY_HASH_SOURCES})
