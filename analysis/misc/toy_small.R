# This is a deprecated toy example, original figure 2
label_mapper = c("2" = "P2", "1" = "P1", "-1" = "T1", "-2" = "T2", "-3" = "T3")
label_mapper_func = function(x) {
  label_mapper[x]
}

extract_mut_counts_history <- function(x_gens, type_graph) {
  if (x_gens$num_gen > 0) {
    m_counts_gens = lapply(1:x_gens$num_gen, function(i) {
      lapply(x_gens$mut_counts_history[[i]], function(x) {
        # singleton vs double
        x[[1]] + x[[2]]
      })
    })
  } else {
    m_counts_gens = list()
  }
  m_counts_gens[[length(m_counts_gens)+1]] = x_gens$end_mut_counts
  names(m_counts_gens) = paste0("type_", x_gens$cell_type, "_gen_", 0:(x_gens$num_gen))
  if (!is.null(x_gens$end_mut_counts_mode)) {
    out_types = type_graph$merge[x_gens$cell_type, ]
    diff_mode_names = paste0("type_", x_gens$cell_type, "_gen_", x_gens$num_gen, "_mode1_lin_", out_types)
    m_counts_gens[[diff_mode_names[1]]] = lapply(x_gens$end_mut_counts_mode, "[[", 1)
    m_counts_gens[[diff_mode_names[2]]] = lapply(x_gens$end_mut_counts_mode, "[[", 2)
  }
  if (!is.null(x_gens$cell_mode_counts)) {
    for (i in 1:x_gens$num_gen) {
      if (x_gens$cell_mode_counts[i+1, 2] > 0) {
        sampled_name = paste0("type_", x_gens$cell_type, "_gen_", i - 1, "_lin_sampled")
        unsampled_name = paste0("type_", x_gens$cell_type, "_gen_", i - 1, "_lin_unsampled")
        unsampled_name_next = paste0("type_", x_gens$cell_type, "_gen_", i, "_unsampled")
        # the sampled has the same muts as the cells for that gen, because unsampled has count zero
        m_counts_gens[[sampled_name]] = lapply(x_gens$mut_counts_history[[i]], function(x) {
          # singleton vs double
          x[[1]] + x[[2]]
        })
        empty_muts = lapply(1:length(x_gens$mut_counts_history[[i]]),
                            function(x) return(matrix(NA, 0, 0)))
        m_counts_gens[[unsampled_name]] = empty_muts
        m_counts_gens[[unsampled_name_next]] = empty_muts
      }
    }

  }
  allele_mats = extract_allele_matrix(m_counts_gens, slot = NULL)
  allele_mats_collapsed = lapply(allele_mats, collapse_reucr_ver)
  allele_mats_collapsed
}

merge_mut_counts <- function(mut_counts_list) {
  # for mergin mut histories from the types
  num_elements = sapply(mut_counts_list, length)
  assertthat::assert_that(is_zero_range(num_elements))
  # type counts
  lapply(1:num_elements[1], function(j) {
    m_alleles = sort(unique(do.call(c, lapply(mut_counts_list, function(x) colnames(x[[j]])))))
    zz = do.call(rbind, lapply(names(mut_counts_list), function(x) {
      m_counts = mut_counts_list[[x]][[j]]
      out = m_counts[, match(m_alleles, colnames(m_counts)), drop = F]
      colnames(out) = m_alleles
      out
    }))
    zz[is.na(zz)] = 0
    zz
  })
}

# template for graph definition
# define nodes and tips
root_id = "2"
node_id  = c("2", "1")
tip_id = c("-1", "-2", "-3")

# define topology
merge0 = do.call(rbind, list("1" = c(-1, -2),
                             "2" = c(1, -3)))
# define unsampled branches (ignoring unsampled for the main figure now)
# trans_sel0 = list("2" = list(list(fraction = c(0.5, 0),
#                                   time = 1.0)))
type_mode_probs0 = list("1" = c(0.4, 0.6, 0, 0, 0),
                        "2" = c(0.75, 0.25, 0, 0, 0))
type_time0 = list("1" = 2.0 + 0.25 * 2,
                  "2" = 2.0)
tip_time = rep(Inf, length(tip_id))
names(tip_time) = tip_id
type_time0 = c(type_time0, tip_time)
type_doubletime = list("1" = 0.25,
                       "2" = 0.4,
                       "-1" = 0.3,
                       "-2" = 0.4,
                       "-3" = 0.2)
toy_graph = make_type_graph(root_id = root_id,
                            node_id = node_id,
                            tip_id = tip_id,
                            merge_matrix = merge0,
                            differntiation_time = type_time0,
                            differntiation_mode_probs = type_mode_probs0,
                            double_time = type_doubletime,
                            founder_size = 1,
                            target_time = 2.1 + 0.25 * 5)
toy_col_mapper <- function(x) {
        col_vec = c("-3" = "#4A6B8A",
                    "1" = "#E699AF",
                    "-1" = "#dd5182",
                    "-2" = "#d73027",
                    "2" = "#D4CB6A")
        return(col_vec[x])
}
g0 = plot_type_graph(toy_graph, toy_col_mapper,
                     node_label_mapper = label_mapper_func,
                     show_node_text = F,
                     effective_time = T)


# aa = sapply(1:800, function(ii) {
set.seed(36)
type_graph = toy_graph
gens0 = make_gens(toy_graph)
gens0 = sample_size_gens(gens0, type_graph, 12)
mut_p = list(mut_rate = 1.,
             recur_prob = 1.,
             recur_vec = c(0.5, 0.2, 0.1, 0.1, 0.05, 0.05),
             active_time = list(c(0, 3.)))
mut_param = do.call(make_mut_param_by_rate, mut_p)
gens0 = simulate_muts(gens0, type_graph, mut_param = mut_param)
mut_counts_type = lapply(gens0, extract_mut_counts_history, type_graph = type_graph)

mut_counts_history = merge_mut_counts(mut_counts_type)
mut_counts_history

g = make_empty_graph()
for(type_gens in gens0[c(tip_id, backward_merge_sequence(type_graph))]) {
        # control for if there is unsampled branches
        if (!is.null(type_gens$cell_mode_counts)) {
                for (i in 0:type_gens$num_gen) {
                        g = g %>% add_vertices(1,
                                               name = paste0("type_", type_gens$cell_type, "_gen_", i),
                                               count = type_gens$cell_mode_counts[i+1, 1],
                                               time = type_gens$start_time + type_gens$double_time * i,
                                               cell_type = type_gens$cell_type,
                                               vertex_type = "cells")
                        if (!is.null(type_gens$sample_size)) {
                                V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[i+1]
                        }
                        if (type_gens$cell_mode_counts[i+1, 2] > 0) {
                                # cell fate vertices
                                g = g %>% add_vertices(1,
                                                       name = paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_sampled"),
                                                       count = type_gens$cell_mode_counts[i+1, 1]/2,
                                                       time = type_gens$start_time + type_gens$double_time * (i - 1),
                                                       cell_type = type_gens$cell_type,
                                                       vertex_type = "Committed cells"
                                ) %>% add_vertices(1,
                                                  name = paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_unsampled"),
                                                  count = type_gens$cell_mode_counts[i+1, 2]/2,
                                                  time = type_gens$start_time + type_gens$double_time * (i - 1),
                                                  cell_type = NA,
                                                  vertex_type = "Committed cells"
                                )
                                if (!is.null(type_gens$sample_size_mode)) {
                                        V(g)[paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_sampled")]$sample_size = type_gens$sample_size_mode[i, 1]
                                        V(g)[paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_unsampled")]$sample_size = 0
                                }
                                g = g %>% add_vertices(1,
                                                       name = paste0("type_", type_gens$cell_type, "_gen_", i, "_unsampled"),
                                                       count = type_gens$cell_mode_counts[i+1, 2],
                                                       time = type_gens$start_time + type_gens$double_time * i,
                                                       cell_type = NA,
                                                       vertex_type = "cells"
                                )
                                V(g)[paste0("type_", type_gens$cell_type, "_gen_", i, "_unsampled")]$sample_size = 0
                        }
                }
                # edges
                for (i in 1:type_gens$num_gen) {
                        if (type_gens$cell_mode_counts[i+1, 2] > 0) {
                                g = g %>%
                                        add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", i - 1),
                                                      paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_unsampled")),
                                                  edge_type = "cell fate commit") %>%
                                        add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", i - 1),
                                                    paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_sampled")),
                                                  edge_type = "cell fate commit") %>%
                                        add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_sampled"),
                                                    paste0("type_", type_gens$cell_type, "_gen_", i)),
                                                  edge_type = "cell division")  %>%
                                        add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", i - 1, "_lin_unsampled"),
                                                    paste0("type_", type_gens$cell_type, "_gen_", i, "_unsampled")),
                                                  edge_type = "cell division")
                        } else {
                                g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", i-1),
                                                      paste0("type_", type_gens$cell_type, "_gen_", i)),
                                                    edge_type = "cell division")
                        }
                }
        } else {
                # when there's no unsampled branching
                # vertices
                for (i in 0:type_gens$num_gen) {
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
                # edges
                for (i in 1:type_gens$num_gen) {
                        g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", i-1),
                                              paste0("type_", type_gens$cell_type, "_gen_", i)),
                                            edge_type = "cell division")
                }
        }
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
# plot(g)

theme0 = theme(
        panel.background = element_rect(fill = NA, color = NA, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=9),
        text = element_text(size = 12)
)
theme1 = theme(
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = NA, size = 10),
        legend.text=element_text(size=9),
        text = element_text(size = 12),
        axis.title.x = element_blank()
)

library(ggraph)
# ggraph(g, layout = "dendrogram", height = time) + geom_node_label(aes(label = name, color = cell_type)) + geom_edge_bend() + ylim(c(type_graph$target_time, 0))
g1 = ggraph(g, layout = "dendrogram", height = time) + geom_edge_bend(aes(linetype = edge_type)) +
  geom_node_point(aes(shape = vertex_type, color = cell_type, size = log2(count+1))) +
        geom_node_text(aes(label = count)) +
        ylim(c(type_graph$target_time, 0)) + scale_size_continuous(range = c(6, 9), guide = F) +
        scale_color_manual(values = toy_col_mapper(sort(unique(V(g)$cell_type))), na.value = "grey") + theme0

V(g)$unsampled_count = V(g)$count - V(g)$sample_size
xy <- ggraph::create_layout(g, "dendrogram", height = time)
V(g)$x <- xy[match(V(g)$name, xy$name), 1]
V(g)$y <- xy[match(V(g)$name, xy$name), 2]

# sample size
# library(scatterpie)
# dd = igraph::as_data_frame(g, "vertices")
# dd$sampled_color = toy_col_mapper(dd$cell_type)
# dd$sampled_color = toy_col_mapper(dd$cell_type)

g2 = ggraph(g, layout = "manual", x = V(g)$x, y = 3 - V(g)$y) +
  geom_edge_bend(aes(linetype = edge_type)) +
  geom_scatterpie(aes(x = x, y =  3 - y, r = log2(count+1)/140+0.05),
                  cols = c("unsampled_count", "sample_size"),
                  data = igraph::as_data_frame(g, "vertices")) + scale_fill_manual(values = c("#d1d1d1", "black")) +
          scale_edge_linetype(guide=F) + theme1


# step 3: mutations
mut_names = c("0", "R-1", "R-2", "R-3", "R-4", "R-5", "R-6")
mut_col_mapper <- function(z) {
        z = stringr::str_replace(z, "\\.", "-")
        c("#D3D3D3", "#98a77b", "#7ba78b", "#7b98a7", "#8b7ba7", "#a77b98", "#a78b7b")[match(z, mut_names)]
}

# single cell
library(data.table)
set.seed(99)

# begin tree with all cells
gens0 = make_gens(toy_graph)
gens0 = sample_size_gens(gens0, type_graph, sample_size = c("-1" = 48, "-2" = 48, "-3" = 32))
# want to try to highlight the sampled cells here, there may be some code from before
# end tree with all cells

gens1 = simulate_sc_muts(gens0,
                         type_graph = type_graph,
                         mut_param = mut_param)
get_gen_name_from_id <- function(id_name) {
        stringr::str_match(id_name, "^(.*?)_(\\d+)$")[, 2]
}

tr = construct_true_lineage(gens1, type_graph$tip_id)
g_tr = as.igraph(tr)
vertices_wide = igraph::as_data_frame(g, what = "vertices")
V(g_tr)$type = get_type_from_id(V(g_tr)$name)
V(g_tr)$time = vertices_wide$time[match(get_gen_name_from_id(V(g_tr)$name), vertices_wide$name)]

# manual mapping time for now
V(g_tr)$time[get_gen_name_from_id(V(g_tr)$name) == "type_1_gen_-1"] = vertices_wide$time[vertices_wide$name == "type_2_gen_4_mode1_lin_1"]
V(g_tr)$time[get_gen_name_from_id(V(g_tr)$name) == "type_-3_gen_-1"] = vertices_wide$time[vertices_wide$name == "type_2_gen_4_mode1_lin_-3"]
V(g_tr)$time[get_gen_name_from_id(V(g_tr)$name) == "type_-1_gen_-1"] = vertices_wide$time[vertices_wide$name == "type_1_gen_1_mode1_lin_-1"]
V(g_tr)$time[get_gen_name_from_id(V(g_tr)$name) == "type_-2_gen_-1"] = vertices_wide$time[vertices_wide$name == "type_1_gen_1_mode1_lin_-2"]

V(g_tr)$time[get_gen_name_from_id(V(g_tr)$name) %in% c("type_-1_gen_2", "type_-2_gen_2", "type_-3_gen_6")] = 3.35

# g_tr = g_tr %>% add_vertices(1, name = "type_2_gen_0_1", type = "2", time = 0) %>% add_edges(c("type_2_gen_0_1",
                                                                                       # "type_2_gen_1_1"))
g4 = ggraph(g_tr, layout = "dendrogram", height = time) + geom_edge_diagonal() +
  geom_node_point(aes(color = type, shape = "single cell"), size = 1.5) +
  scale_size_continuous(guide = F) +
        scale_shape_manual(values = 15) +
        scale_color_manual(values = toy_col_mapper(unique(V(g_tr)$type)),
                           guide = F) + ylim(c(type_graph$target_time, 0)) + theme0

sc_muts = get_barcodes(gens1, type_graph$tip_id)
sc_muts = remove_allele_ver(sc_muts)
all(c("0", "R-1", "R-2", "R-3", "R-4", "R-5", "R-6") %in% sc_muts[, 1])


get_barcodes_history <- function(gens1, type_graph) {
        do.call(rbind, lapply(c(forward_merge_sequence(type_graph), type_graph$tip_id), function(x) {
                do.call(rbind, c(lapply(gens1[[x]]$barcode_history,
                                        function(z) rbind(z[[1]], z[[2]])),
                                 list(gens1[[x]]$end_barcodes)))
        }))
}
get_commit_barcodes <- function(gens1, type_graph) {
        do.call(rbind, lapply(forward_merge_sequence(type_graph), function(x) {
                x_gens = gens1[[x]]
                out = NULL
                if (x_gens$active) {
                  out_types = type_graph$merge[x_gens$cell_type, ]
                  b_dist = x_gens$end_barcodes_mode
                  out = rbind(out, do.call(rbind, lapply(1:2, function(j) {
                    b_mat = b_dist[[j]]
                    name_split = stringr::str_match(rownames(b_mat), "^(type_.*?_gen_\\d+)_(.*)$")
                    rownames(b_mat) = paste0(name_split[, 2], "_mode1_lin_", out_types[j], "_", name_split[, 3])
                    b_mat
                  })))
                }
                if (!is.null(x_gens$cell_mode_counts)) {
                        trans_gens = which(x_gens$cell_mode_counts[, 2] > 0) - 2
                        out = rbind(out, do.call(rbind, lapply(trans_gens, function(i) {
                                b_dist = x_gens$barcode_history[[i+1]]
                                b_mat = rbind(b_dist[[1]], b_dist[[2]])
                                name_split = stringr::str_match(rownames(b_mat), "^(type_.*?_gen_\\d+)_(.*)$")
                                rownames(b_mat) = paste0(name_split[, 2], "_lin_sampled", "_", name_split[, 3])
                                b_mat
                        })))
                }
                out
        }))
}

barcodes_his = remove_allele_ver(get_barcodes_history(gens1, type_graph))
barcodes_commit = remove_allele_ver(get_commit_barcodes(gens1, type_graph))
# type_gen_cell_total = table(c(get_gen_name_from_id(rownames(barcodes_his)),
#                                  get_gen_name_from_id(rownames(barcodes_commit))))
# all(type_gen_cell_total == vertices_wide$sample_size[match(names(type_gen_cell_total), vertices_wide$name)])

all(V(g_tr)$name %in% rownames(barcodes_his))
V(g_tr)$mut = barcodes_his[match(V(g_tr)$name, rownames(barcodes_his)), 1]
g5 = ggraph(g_tr, layout = "dendrogram", height = time) +geom_edge_diagonal() + geom_node_point(aes(color = mut, shape = "single cell"), size = 1.5) +
  # geom_node_point(aes(color = mut, shape = "single cell"), size = 3) +
        scale_size_continuous(guide = F) + scale_color_manual(values = mut_col_mapper(sort(unique(V(g_tr)$mut)))) +
        scale_shape_manual(values = 15) +
  ylim(c(type_graph$target_time, 0)) + theme0

# library(tidyverse)
# vertices_wide <- igraph::as_data_frame(g, "vertices")
# mut_counts_df = data.frame(mut_counts_history[[1]])
# mut_counts_df$name = rownames(mut_counts_df)
# vertices_wide <- merge(vertices_wide, mut_counts_df, by = "name")

# sum up single cell for bulk example
barcodes_pooled = do.call(rbind, c(lapply(split(barcodes_his, get_gen_name_from_id(rownames(barcodes_his))), function(x) table(x)[mut_names]),
  lapply(split(barcodes_commit, get_gen_name_from_id(rownames(barcodes_commit))), function(x) table(x)[mut_names])))
barcodes_pooled[is.na(barcodes_pooled)] = 0
colnames(barcodes_pooled) = mut_names
barcodes_pooled = as.data.frame(barcodes_pooled)
barcodes_pooled$name = rownames(barcodes_pooled)
vertices_wide <- left_join(vertices_wide, barcodes_pooled, by = "name")
vertices_wide[is.na(vertices_wide)] = 0

vertices_long <- vertices_wide[, colnames(vertices_wide) %in% c("name", mut_names)] %>%
        gather("mut", "value", -name)
# vertices_long$mut[vertices_long$mut == "X0"] = "0"
# create the bar charts
bar_list <- lapply(V(g)$name, function(x) {
        gt_plot <- ggplotGrob(
                ggplot(vertices_long[vertices_long$name == x, ]) +
                        geom_col(
                                aes(factor(name), value, group = name, fill = mut),

                                color = NA
                        ) +
                        labs(x = NULL, y = NULL) +
                        coord_cartesian(expand = FALSE) + coord_flip() +
                        scale_fill_manual(values = mut_col_mapper(unique(vertices_long$mut))) +
                        theme(
                                legend.position = "none",
                                panel.background = element_rect(fill = "white", colour = NA),
                                line = element_blank(),
                                text = element_blank()
                        )
        )
        panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel", ]
        gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
})
names(bar_list) = V(g)$name

# convert to custom annotation
annot_list <- lapply(V(g)$name, function(x) {
        vertex_attr = vertices_wide[vertices_wide$name == x, ]
        xmin <- vertex_attr$x - (log2(vertex_attr$sample_size+1)/30 + 0.015)
        xmax <- vertex_attr$x + (log2(vertex_attr$sample_size+1)/30 + 0.015)
        ymin <- vertex_attr$y - (log2(vertex_attr$sample_size+1)/30 + 0.015)
        ymax <- vertex_attr$y + (log2(vertex_attr$sample_size+1)/30 + 0.015)
        annotation_custom(
                bar_list[[x]],
                xmin = xmin,
                xmax = xmax,
                ymin = 3 - ymin,
                ymax = 3 - ymax
        )
})
p = ggraph(g, layout = "manual", x = V(g)$x, y = 3 - V(g)$y) + scale_edge_linetype(guide = F) +
        geom_edge_bend(aes(linetype = edge_type)) + theme1
g3 = Reduce("+", annot_list, p)
library(patchwork)

(g0 + g1 + g2 + g4) + plot_layout(nrow = 1, widths = c(1, 2, 2, 3), guides = "collect")


g_all = (g0 | g1) / ( g2 | g4) /
  (g3 | g5) + plot_layout(guides = "collect")
push_pdf(g_all, "toy_all_v1", dir = "./plots/toy/", width = 7.2, height = 7.2)

