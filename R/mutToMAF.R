#' NGS mutation data perform
#'
#' This function converts mutation data from NGS sequencing into MAF format files and performs filtering.
#'
#' @param root_dir            生信拉的数据存储的文件路径.
#' @param panel               测序的panel，例如639
#' @param clin                data.frame格式入组患者编号，表头名称需为Tumor_Sample_Barcode
#' @param saveDATA            if TRUE 将在root_dir路径下新建一个SNV文件夹并生成文件，默认FALSE
#' @param mut_filter          突变筛选。defalut mut_filter=TRUE
#' @param tumor_t             该位点肿瘤样本reads数
#' @param site_depth          该位点测序深度
#' @param hotspot_vaf         热点突变丰度.Druggable为'0|1|2|-'标签
#' @param non_hotspot_vaf     非热点突变丰度
#' @param hotspotloss_vaf     热点失活突变丰度
#' @param non_hotspotloss_vaf 非热点失活突变丰度
#' @param loss_genes          抑癌基因/HRR/MMR相关基因，和失活突变筛选有关
#' @param tumor               对应的癌种tumor：运行data("hotsgenes"),colnames(hotsgenes)
#' @param hots_genes          癌种对应的用药等级1，2和R1的基因。例如肺癌：hots_genes=c("ALK","BRAF","EGFR","ERBB2","KRAS","MET","RET","ROS1")
#' @param artifact            假突变标签，默认artifact=c("Germline","background","base_qual","normal_artifact","weak_evidence","contamination","strand_bias","haplotype","multiallelic","chemotherapy_site","suspected_bg","low_allele_frac")
#' @return A MAF type data for NGS mutation data with filtering
#' @export
##' @author Pengmin Yang
#' @examples
#' data("clindata")
#' root_dir <-system.file("example",  package = "Mypackage")
#' result <-mutToMAF(root_dir=root_dir,clin=clindata,saveDATA=FALSE,mut_filter=TRUE,
#'                   tumor_t=10,site_depth=100,hotspot_vaf=0.009,
#'                  non_hotspot_vaf=0.045,hotspotloss_vaf=0.095,non_hotspotloss_vaf=0.195)
#' print(result)
mutToMAF<-function (root_dir, clin, saveDATA = FALSE, mut_filter = TRUE,
                    tumor_t = 10, site_depth = 100, hotspot_vaf = 0.009, non_hotspot_vaf = 0.045,
                    hotspotloss_vaf = 0.095, non_hotspotloss_vaf = 0.195, panel = 639,
                    loss_genes = c("APC", "ARID1A", "ARID2", "ATM", "ATR", "ATRX",
                                   "BAP1", "BARD1", "BLM", "BRCA1", "BRCA2", "BRIP1", "CDK12",
                                   "CDKN2A", "CHEK1", "CHEK2", "CREBBP", "DNMT3A", "FAM175A",
                                   "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCL",
                                   "FANCM", "FBXW7", "FLCN", "GEN1", "GNA13", "KDM5C", "KDM6A",
                                   "KMT2A", "KMT2C", "KMT2D", "MAP3K1", "MLH1", "MRE11",
                                   "MSH2", "MSH3", "MSH6", "MYC", "NBN", "NF1", "NF2", "NOTCH1",
                                   "NSD1", "PALB2", "PBRM1", "PMS2", "POLD1", "PPP2R1A",
                                   "PPP2R2A", "PRKDC", "PTCH1", "PTEN", "RAD50", "RAD51",
                                   "RAD51B", "RAD51C", "RAD51D", "RAD54L", "RB1", "RECQL",
                                   "RINT1", "RNF43", "SDHA", "SDHB", "SDHC", "SDHD", "SETD2",
                                   "SLX4", "SMAD4", "SMARCA4", "SMARCB1", "STK11", "TP53",
                                   "TSC1", "TSC2", "VHL", "WRN", "XRCC2"), tumor = NULL,
                    hots_genes = c("ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET",
                                   "RET", "ROS1"), artifact = c("Germline", "background",
                                                                "base_qual", "normal_artifact", "weak_evidence", "contamination",
                                                                "strand_bias", "haplotype", "multiallelic", "chemotherapy_site",
                                                                "suspected_bg", "low_allele_frac"))
{
  library(dplyr)
  library(stringr)
  library(Hmisc)
  merged_df <- data.frame()
  panel <- paste0("_", panel, "_")
  files <- list.files(path = root_dir, full.names = TRUE)
  count <- 0
  start_time <- Sys.time()
  for (file_path in files) {
    loop_start_time <- Sys.time()
    if (grepl("gzy", file_path) && grepl(paste(panel, collapse = "|"),
                                         file_path) && !grepl("dedup", file_path) && grepl("review_for_report\\.xls$",
                                                                                           file_path) && !grepl("X_g|PCR", file_path)) {
      df <- read.table(file_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, colClasses = c(CHROM = "character",
                                                                SF_TVAF = "numeric", Druggable = "character"),
                       na.strings = "*")
      df <- df[, c(1, 3:4, 6, 10:17, 19:20)]
      df$DP <- NA
      df$T1 <- NA
      df$depth <- NA
      for (j in 1:nrow(df)) {
        split_df <- strsplit(df[j, "DP_T"], split = "|",
                             fixed = TRUE)
        df[j, "DP"] <- as.numeric(split_df[[1]][[1]])
        df[j, "T1"] <- as.numeric(split_df[[1]][[2]])
        df[j, "depth"] <- as.numeric(split_df[[1]][[1]]) +
          as.numeric(split_df[[1]][[2]])
      }
      df <- df[, -4]
      if (nrow(df) > 0) {
        file_name <- basename(file_path)
        df$file_name <- file_name
        if (grepl("ct", file_name)) {
          Tumor_Sample_Barcode <- gsub("gzy|_ct.*", "",
                                       file_name)
          Tumor_Sample_Barcode <- sub("([A-Z].*)", "",
                                      Tumor_Sample_Barcode)
          df$Tumor_Sample_Barcode <- Tumor_Sample_Barcode
        }
        else if (grepl("Q", file_name) && !grepl("C",
                                                 file_name)) {
          Tumor_Sample_Barcode <- gsub("gzy|Q.*", "",
                                       file_name)
          Tumor_Sample_Barcode <- sub("([A-Z].*)", "",
                                      Tumor_Sample_Barcode)
          df$Tumor_Sample_Barcode <- Tumor_Sample_Barcode
        }
        else if (grepl("C", file_name)) {
          Tumor_Sample_Barcode <- gsub("gzy|C.*", "",
                                       file_name)
          df$Tumor_Sample_Barcode <- Tumor_Sample_Barcode
        }
        else {
          Tumor_Sample_Barcode <- gsub("gzy|_g.*", "",
                                       file_name)
          Tumor_Sample_Barcode <- sub("([A-Z].*)", "",
                                      Tumor_Sample_Barcode)
          df$Tumor_Sample_Barcode <- Tumor_Sample_Barcode
        }

        if (Tumor_Sample_Barcode %in% clin[, 1]) {
          colnames(df)[colnames(df) == "Gene"] <- "Hugo_Symbol"
          df$pHGVS <- ifelse(df$pHGVS == "", NA, df$pHGVS)
          if (mut_filter) {
            df$T1 <- as.numeric(df$T1)
            df$Druggable <- as.character(df$Druggable)
            df$Tumor_VAF <- as.numeric(df$Tumor_VAF)
            df$pHGVS <- as.character(df$pHGVS)
            df$FILTER <- as.character(df$FILTER)
            df$depth <- as.numeric(df$depth)
            if (is.null(hots_genes)) {
              data("hotsgenes",package="Mypackage")
              if (length(tumor) == 1) {
                if (!is.null(hotsgenes[[tumor]]) && is.character(hotsgenes[[tumor]])) {
                hots_genes <- strsplit(hotsgenes[[tumor]],
                                       split = "|", fixed = TRUE)[[1]]
                }else {
                  warning(paste(tumor,"is NULL or not a character vector."))
                }
              }
              else if (length(tumor) > 1) {
                hots_genes <- unlist(lapply(tumor, function(t) {
                  if (!is.null(hotsgenes[[t]]) && is.character(hotsgenes[[t]])) {
                  strsplit(hotsgenes[[t]], split = "|",
                           fixed = TRUE)[[1]]
                  }else {
                    warning(paste(t,"is NULL or not a character vector."))
                  }

                }))
                hots_genes <- unique(hots_genes)
              }
            }
            else {
              hots_genes <- hots_genes
            }
            loss_genes <- loss_genes
            df <- df %>% filter(T1 >= tumor_t, depth >=
                                  site_depth) %>% mutate(is_hotspot = grepl("0|1|2|-",
                                                                            Druggable), is_hotsgene = Hugo_Symbol %in%
                                                           hots_genes, is_loss_gene = Hugo_Symbol %in%
                                                           loss_genes, is_truncating = grepl("fs|\\*",
                                                                                             pHGVS), artifact = grepl(paste(artifact,
                                                                                                                            collapse = "|"), FILTER)) %>% mutate(mutation_class = case_when(is_hotspot &
                                                                                                                                                                                              is_hotsgene & !is_truncating & Tumor_VAF >=
                                                                                                                                                                                              hotspot_vaf ~ 1, !artifact & is_hotspot &
                                                                                                                                                                                              !is_hotsgene & !is_truncating & Tumor_VAF >=
                                                                                                                                                                                              hotspot_vaf ~ 1.1, !artifact & is_hotspot &
                                                                                                                                                                                              is_loss_gene & is_truncating & Tumor_VAF >=
                                                                                                                                                                                              hotspotloss_vaf ~ 2, !artifact & is_hotspot &
                                                                                                                                                                                              !is_loss_gene & is_truncating & Tumor_VAF >=
                                                                                                                                                                                              non_hotspotloss_vaf ~ 3, !artifact & !is_hotspot &
                                                                                                                                                                                              !is_truncating & Tumor_VAF >= non_hotspot_vaf &
                                                                                                                                                                                              Tumor_VAF <= 0.9 ~ 4, !artifact & !is_hotspot &
                                                                                                                                                                                              is_truncating & Tumor_VAF >= non_hotspotloss_vaf &
                                                                                                                                                                                              Tumor_VAF <= 0.9 ~ 5)) %>% filter(!is.na(mutation_class)) %>%
              select(-c(is_hotspot, is_loss_gene, is_truncating,
                        artifact))
            df <- unique(df)
            df <- df[!(df$Hugo_Symbol == "PTEN" & substr(df$cHGVS,
                                                         1, 5) == "c.802"), ]
          }
          if (nrow(df) > 0) {
            colnames(df)[colnames(df) == "CHROM"] <- "Chromosome"
            colnames(df)[colnames(df) == "Start"] <- "Start_Position"
            colnames(df)[colnames(df) == "End"] <- "End_Position"
            colnames(df)[colnames(df) == "REF"] <- "Reference_Allele"
            colnames(df)[colnames(df) == "ALT"] <- "Tumor_Seq_Allele2"
            colnames(df)[colnames(df) == "Transcript"] <- "RefSeq"
            colnames(df)[colnames(df) == "DP"] <- "t_ref_count"
            colnames(df)[colnames(df) == "T1"] <- "t_alt_count"
            split_columns <- strsplit(df$Report_C2, ",")
            num_columns <- sapply(split_columns, length)
            last_values <- sapply(split_columns, function(x) ifelse(length(x) >
                                                                      0, x[length(x)], ""))
            last_values[!grepl("p\\.", last_values)] <- NA
            df$HGVSp_Short <- last_values
            df$Variant_Type <- ifelse(nchar(df$Tumor_Seq_Allele2) ==
                                        nchar(df$Reference_Allele), "SNP", ifelse(nchar(df$Tumor_Seq_Allele2) >
                                                                                    nchar(df$Reference_Allele), "INS", "DEL"))
            df$Variant_Classification <- NA
            for (i in 1:nrow(df)) {
              if (df[i, "Variant_Type"] == "SNP") {
                if (!is.na(df[i, "pHGVS"])) {
                  if ((((substr(df[i, "pHGVS"], nchar(df[i,
                                                         "pHGVS"]), nchar(df[i, "pHGVS"])) !=
                         "*") || (grepl(">", df[i, "pHGVS"]) &&
                                  !grepl("\\*", df[i, "pHGVS"]))) &&
                       (substr(df[i, "pHGVS"], nchar(df[i,
                                                        "pHGVS"]), nchar(df[i, "pHGVS"])) !=
                        "=") && (substr(df[i, "pHGVS"],
                                        3, 3) != "*") && (substr(df[i, "pHGVS"],
                                                                 3, 5) != "Met"))) {
                    df[i, "Variant_Classification"] <- "Missense_Mutation"
                  }
                  else if ((((substr(df[i, "pHGVS"],
                                     nchar(df[i, "pHGVS"]), nchar(df[i,
                                                                     "pHGVS"])) == "*") || (grepl(">",
                                                                                                  df[i, "pHGVS"]) && grepl("\\*", df[i,
                                                                                                                                     "pHGVS"]))) && (substr(df[i, "pHGVS"],
                                                                                                                                                            nchar(df[i, "pHGVS"]), nchar(df[i,
                                                                                                                                                                                            "pHGVS"])) != "=") && (substr(df[i,
                                                                                                                                                                                                                             "pHGVS"], 3, 3) != "*") && (substr(df[i,
                                                                                                                                                                                                                                                                   "pHGVS"], 3, 5) != "Met"))) {
                    df[i, "Variant_Classification"] <- "Nonsense_Mutation"
                  }
                  else if ((substr(df[i, "pHGVS"], nchar(df[i,
                                                            "pHGVS"]), nchar(df[i, "pHGVS"])) ==
                            "=") && (substr(df[i, "pHGVS"], 3,
                                            3) != "*") && (substr(df[i, "pHGVS"],
                                                                  3, 5) != "Met")) {
                    df[i, "Variant_Classification"] <- "Silent"
                  }
                  else if ((substr(df[i, "pHGVS"], 3,
                                   3) == "*") && (substr(df[i, "pHGVS"],
                                                         3, 5) != "Met")) {
                    df[i, "Variant_Classification"] <- "Nonstop_Mutation"
                  }
                  else if (substr(df[i, "pHGVS"], 3,
                                  5) == "Met") {
                    df[i, "Variant_Classification"] <- "Translation_Start_Site"
                  }
                }
                else if (is.na(df[i, "pHGVS"])) {
                  if ((grepl(">", df[i, "cHGVS"])) &&
                      (!grepl("del|ins", df[i, "cHGVS"])) &&
                      (!grepl("_", df[i, "cHGVS"]))) {
                    num <- str_extract(df[i, "cHGVS"],
                                       "(?<=\\+|-)(\\d+)(?=[A-Za-z])")
                    num <- ifelse(is.na(num), 0, num)
                    num <- as.numeric(num)
                    if (!is.na(num)) {
                      if (num <= 2) {
                        df[i, "Variant_Classification"] <- "Splice_Site"
                      }
                      else if (num > 2 && num <= 10) {
                        df[i, "Variant_Classification"] <- "Splice_Region"
                      }
                      else if (num > 10) {
                        df[i, "Variant_Classification"] <- "Intron"
                      }
                    }
                  }
                  else if ((!grepl(">", df[i, "cHGVS"])) &&
                           (grepl("del|ins", df[i, "cHGVS"]) &&
                            (grepl("_", df[i, "cHGVS"])))) {
                    num1 <- str_extract(df[i, "cHGVS"],
                                        "(?<=\\+|-)(\\d+)(?=_)")
                    num1 <- ifelse(is.na(num1), 0, num1)
                    num2 <- str_extract(df[i, "cHGVS"],
                                        "(?<=\\+|-)(\\d+)(?=del)")
                    num2 <- ifelse(is.na(num2), 0, num2)
                    num1 <- as.numeric(num1)
                    num2 <- as.numeric(num2)
                    if ((!is.na(num1)) && ((!is.na(num2)))) {
                      if ((num1 <= 2 || num2 <= 2) ||
                          (num1 <= 2 && num2 <= 2)) {
                        df[i, "Variant_Classification"] <- "Splice_Site"
                      }
                      else if (((num1 > 2 && num1 < 10) ||
                                (num2 > 2 && num2 < 10)) || ((num1 >
                                                              2 && 10 < num1) && (num2 > 2 &&
                                                                                  num2 < 10))) {
                        df[i, "Variant_Classification"] <- "Splice_Region"
                      }
                      else if (num1 > 10 && num2 > 10) {
                        df[i, "Variant_Classification"] <- "Intron"
                      }
                    }
                  }
                }
              }
              else if (df[i, "Variant_Type"] == "DEL") {
                if (!is.na(df[i, "pHGVS"]) && grepl("fs",
                                                    df[i, "pHGVS"])) {
                  df[i, "Variant_Classification"] <- "Frame_Shift_Del"
                }
                else if (!is.na(df[i, "pHGVS"]) && (grepl("\\*$",
                                                          df[i, "pHGVS"]))) {
                  df[i, "Variant_Classification"] <- "Frame_Shift_Del"
                }
                else if (!is.na(df[i, "pHGVS"]) && grepl(">",
                                                         df[i, "pHGVS"]) && grepl(">[^>]*\\*$",
                                                                                  df[i, "pHGVS"])) {
                  df[i, "Variant_Classification"] <- "Frame_Shift_Del"
                }
                else if ((is.na(df[i, "pHGVS"])) && (!is.na(df[i,
                                                               "cHGVS"])) && (!grepl("_", df[i, "cHGVS"])) &&
                         (grepl("del", df[i, "cHGVS"]))) {
                  num3 <- str_extract(df[i, "cHGVS"],
                                      "(?<=\\+|-)(\\d+)(?=del)")
                  num3 <- ifelse(is.na(num3), 0, num3)
                  num3 <- as.numeric(num3)
                  if (num3 <= 2) {
                    df[i, "Variant_Classification"] <- "Frame_Shift_Del"
                  }
                  else if (num3 > 2 && num3 <= 10) {
                    df[i, "Variant_Classification"] <- "Splice_Region"
                  }
                  else if (num3 > 10) {
                    df[i, "Variant_Classification"] <- "Intron"
                  }
                }
                else if (!is.na(df[i, "pHGVS"]) && (grepl("del|ins",
                                                          df[i, "pHGVS"])) && grepl("_", df[i,
                                                                                            "cHGVS"])) {
                  df[i, "Variant_Classification"] <- "In_Frame_Del"
                }
                else if ((is.na(df[i, "pHGVS"])) && (!is.na(df[i,
                                                               "cHGVS"])) && (grepl("_", df[i, "cHGVS"])) &&
                         (grepl("del", df[i, "cHGVS"]))) {
                  num4 <- str_extract(df[i, "cHGVS"],
                                      "(?<=\\+|-)(\\d+)(?=_)")
                  num4 <- ifelse(is.na(num4), 0, num4)
                  num5 <- str_extract(df[i, "cHGVS"],
                                      "(?<=\\+|-)(\\d+)(?=del)")
                  num5 <- ifelse(is.na(num5), 0, num5)
                  num4 <- as.numeric(num4)
                  num5 <- as.numeric(num5)
                  if (((num4 >= 3 && num4 <= 10) || (num5 >=
                                                     3 && num5 <= 10)) || ((num4 >= 3 &&
                                                                            num4 <= 10) && (num5 >= 3 && num5 <=
                                                                                            10))) {
                    df[i, "Variant_Classification"] <- "Splice_Region"
                  }
                  else if (num4 > 10 && num5 > 10) {
                    df[i, "Variant_Classification"] <- "Intron"
                  }
                  else if (((num4 <= 2 || num5 <= 2) ||
                            (num4 <= 2 && num5 <= 2))) {
                    df[i, "Variant_Classification"] <- "Frame_Shift_Del"
                  }
                }
                else if (is.na(df[i, "pHGVS"]) && !is.na(df[i,
                                                            "cHGVS"]) && grepl("\\*", df[i, "cHGVS"])) {
                  df[i, "Variant_Classification"] <- "Nonstop_Mutation"
                }
                else {
                  df[i, "Variant_Classification"] <- "Unknown_del"
                }
              }
              else if (df[i, "Variant_Type"] == "INS") {
                if (!is.na(df[i, "pHGVS"]) && grepl("fs",
                                                    df[i, "pHGVS"])) {
                  df[i, "Variant_Classification"] <- "Frame_Shift_Ins"
                }
                else if (!is.na(df[i, "pHGVS"]) && grepl("\\*$",
                                                         df[i, "pHGVS"])) {
                  df[i, "Variant_Classification"] <- "Frame_Shift_Ins"
                }
                else if (!is.na(df[i, "pHGVS"]) && grepl(">",
                                                         df[i, "pHGVS"]) && grepl(">[^>]*\\*$",
                                                                                  df[i, "pHGVS"])) {
                  df[i, "Variant_Classification"] <- "Frame_Shift_Ins"
                }
                else if (!is.na(df[i, "pHGVS"]) && grepl("ins|dup",
                                                         df[i, "pHGVS"])) {
                  df[i, "Variant_Classification"] <- "In_Frame_Ins"
                }
                else if (is.na(df[i, "pHGVS"]) && (grepl("del|dup|ins",
                                                         df[i, "cHGVS"])) && (!grepl("_", df[i,
                                                                                             "cHGVS"]))) {
                  num6 <- str_extract(df[i, "cHGVS"],
                                      "(?<=\\+|-)(\\d+)(?=del|ins|dup)")
                  num6 <- ifelse(is.na(num6), 0, num6)
                  num6 <- as.numeric(num6)
                  if (num6 <= 2) {
                    df[i, "Variant_Classification"] <- "Frame_Shift_Ins"
                  }
                  else if (num6 > 2 && num6 <= 10) {
                    df[i, "Variant_Classification"] <- "Splice_Region"
                  }
                  else if (num6 > 10) {
                    df[i, "Variant_Classification"] <- "Intron"
                  }
                }
                else if (is.na(df[i, "pHGVS"]) && (grepl("del|dup|ins",
                                                         df[i, "cHGVS"])) && (grepl("_", df[i,
                                                                                            "cHGVS"]))) {
                  num7 <- str_extract(df[i, "cHGVS"],
                                      "(?<=\\+|-)(\\d+)(?=_)")
                  num7 <- ifelse(is.na(num7), 0, num7)
                  num8 <- str_extract(df[i, "cHGVS"],
                                      "(?<=\\+|-)(\\d+)(?=del|ins|dup)")
                  num8 <- ifelse(is.na(num8), 0, num8)
                  num7 <- as.numeric(num7)
                  num8 <- as.numeric(num8)
                  if ((num7 >= 3 && num7 <= 10) && (num8 >=
                                                    3 && num8 <= 10)) {
                    df[i, "Variant_Classification"] <- "Splice_Region"
                  }
                  else if ((num7 >= 3 && num8 > 10) ||
                           (num7 > 10 && num8 >= 3)) {
                    df[i, "Variant_Classification"] <- "Intron"
                  }
                  else if (num7 > 10 && num8 > 10) {
                    df[i, "Variant_Classification"] <- "Intron"
                  }
                  else if (((num7 <= 2 || num8 <= 2) ||
                            (num7 <= 2 && num8 <= 2))) {
                    df[i, "Variant_Classification"] <- "Frame_Shift_Ins"
                  }
                }
                else {
                  df[i, "Variant_Classification"] <- "Unknown_Ins"
                }
              }
              else {
                df[i, "Variant_Classification"] <- "Unknown"
              }
            }
          }
          colnames(df)[colnames(df) == "cHGVS"] <- "HGVSc"
          colnames(df)[colnames(df) == "pHGVS"] <- "HGVSp"
          if (nrow(df) > 0) {
            merged_df <- bind_rows(merged_df, df)
          }
          else {
            merged_df <- merged_df
          }
        }
      }
    }
    loop_end_time <- Sys.time()
    loop_time <- loop_end_time - loop_start_time
    if (Tumor_Sample_Barcode %in% clin[, 1]) {
      cat("\033[31m", "Loop", file_name, "time:", paste(loop_time,
                                                        "s"), "\n", "\033[0m\n")
    }
    else {
      cat("Loop", file_name, "time:", paste(loop_time,
                                            "s"), "\n")
    }
    count <- count + 1
  }
  end_time <- Sys.time()
  total_time <- end_time - start_time
  cat("Total:", count, "\n")
  if (saveDATA) {
    output_dir <- file.path(root_dir, "SNV")
    if (!file.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_file <- file.path(output_dir, "SNV.csv")
    write.csv(merged_df, output_file, row.names = FALSE,
              fileEncoding = "UTF-8")
  }
  labels_list <- list(
    "Tumor_Sample_Barcode" = "患者ID",
    "file_name" = "患者NGS检测输出文件",
    "Variant_Type" = "突变类型",
    "Variant_Classification" = "突变亚类",
    "Druggable" ="用药等级",
    "Tumor_VAF"="突变丰度",
    "Chromosome"="染色体",
    "Start_Position"="突变起始位点",
    "End_Position"="突变终末位点",
    "Reference_Allele"="参考碱基",
    "Tumor_Seq_Allele2" ="突变碱基",
    "FILTER"="突变clinvar标签",
    "Hugo_Symbol"="基因symbol",
    "RefSeq"="转录本",
    "HGVSc"="碱基突变",
    "HGVSp"="氨基酸突变",
    "t_ref_count"="未突变reads数",
    "t_alt_count" ="突变reads数",
    "depth"="测序深度",
    "HGVSp_Short"="氨基酸突变（短码）",
    "is_hotsgene" ="该癌种对应的level 1，2和R1用药等级的基因",
    "mutation_class"="突变等级"
  )
  for (col in names(labels_list)) {
    label(merged_df[[col]]) <- labels_list[[col]]
  }
  return(merged_df)
}
