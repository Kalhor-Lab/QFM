read_id_counts <- function(fn) {
        out = read_tsv(fn, col_names = F)
        if (nrow(out) > 0) {
                out = out %>% rename(id = X1, count = X2)
                out = select(out, id, count)
                return(out)
        } else {
                return(NULL)
        }
}

read_id_spacer_counts <- function(fn) {
        if (!file.exists(fn)) {
                return(NULL)
        }
        out = read_tsv(fn, col_names = F)
        if (nrow(out) > 0) {
                out = out %>% rename(id = X1, spacer = X2, count = X3)
                out = select(out, id, spacer, count)
                return(out)
        } else {
                return(NULL)
        }
}
merge_id_spacer_counts <- function(x, y) {
        if (is.null(x) & !is.null(y)) {
                return(y)
        }
        if (is.null(y) & !is.null(x)) {
                return(x)
        }
        if (is.null(x) & is.null(y)) {
                return(NULL)
        }
        full_join(x, y, by = c("id", "spacer")) %>% mutate(count.x = replace_na(count.x, 0),
                                                           count.y = replace_na(count.y, 0)) %>%
                mutate(count = count.x + count.y) %>%
                select(id, spacer, count) %>% arrange(desc(count))
}

merge_id_counts <- function(x, y) {
        if (is.null(x) & !is.null(y)) {
                return(y)
        }
        if (is.null(y) & !is.null(x)) {
                return(x)
        }
        if (is.null(x) & is.null(y)) {
                return(NULL)
        }
        full_join(x, y, by = "id") %>% mutate(count.x = replace_na(count.x, 0),
                                              count.y = replace_na(count.y, 0)) %>%
                mutate(count = count.x + count.y) %>%
                select(id, count) %>% arrange(desc(count))
}

merge_spacer_counts <- function(x, y) {
        if (is.null(x) & !is.null(y)) {
                return(y)
        }
        if (is.null(y) & !is.null(x)) {
                return(x)
        }
        if (is.null(x) & is.null(y)) {
                return(NULL)
        }
        full_join(x, y, by = c("spacer")) %>% mutate(count.x = replace_na(count.x, 0),
                                                     count.y = replace_na(count.y, 0)) %>%
                mutate(count = count.x + count.y) %>%
                select(-c(count.x, count.y)) %>%
                arrange(desc(count))
}
