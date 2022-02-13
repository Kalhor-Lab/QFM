# new toy figure that doesn't include intermediate generations

j = 9
type_graph = all_graphs[[exp_params$big_graph_id[j]]]
tip_id = type_graph$tip_id

g = make_empty_graph()
for(type_gens in gens0[c(tip_id, backward_merge_sequence(type_graph))]) {
        message(type_gens$cell_type)
        # ignore cases of unsampled branches for now
        # when there's no unsampled branching
        # vertices
        for (i in c(0, type_gens$num_gen)) {
                g = g %>% add_vertices(1,
                                       name = paste0("type_", type_gens$cell_type, "_gen_", i),
                                       count = type_gens$cell_counts[i+1],
                                       time = type_gens$start_time + type_gens$double_time * i,
                                       cell_type =type_gens$cell_type,
                                       vertex_type = "cells")
                if (!is.null(type_gens$sample_size)) {
                        V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[i+1]
                }
        }
        # one edge connecting first and last gen
        g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", 0),
                              paste0("type_", type_gens$cell_type, "_gen_", type_gens$num_gen)),
                            edge_type = "cell division")
        # cell type transition vertives and edges
        if (type_gens$active) {
                out_types = type_graph$merge[type_gens$cell_type, ]
                diff_mode_names = paste0("type_", type_gens$cell_type, "_gen_", type_gens$num_gen, "_mode1_lin_", out_types)
                g =  g %>% add_vertices(1,
                                        name = diff_mode_names[1],
                                        count = type_gens$end_mode_counts[1],
                                        time = type_gens$end_time,
                                        cell_type =type_gens$cell_type,
                                        vertex_type = "Committed cells") %>%
                        add_vertices(1,
                                     name = diff_mode_names[2],
                                     count = type_gens$end_mode_counts[2],
                                     time = type_gens$end_time,
                                     cell_type =type_gens$cell_type,
                                     vertex_type = "Committed cells")
                if (!is.null(type_gens$end_mode_sample_size)) {
                        V(g)[diff_mode_names[1]]$sample_size = type_gens$end_mode_sample_size[1]
                        V(g)[diff_mode_names[2]]$sample_size = type_gens$end_mode_sample_size[2]
                }
                g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", type_gens$num_gen),
                                      diff_mode_names[1]),
                                    edge_type = "cell fate commit") %>%
                        add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", type_gens$num_gen),
                                    diff_mode_names[2]),
                                  edge_type = "cell fate commit")
                out_gen_names = paste0("type_", out_types, "_gen_0")
                g = g %>% add_edges(c(diff_mode_names[1],
                                      out_gen_names[1]),
                                    edge_type = "cell division") %>%
                        add_edges(c(diff_mode_names[2],
                                    out_gen_names[2]),
                                  edge_type = "cell division")
        }
}
