#' @title Aggregate the PSM level results from multiple search algorithms
#'
#' @param saveas a character string specifying the file names that scores the results of APIR.
#' APIR outputs three files: a PSM level output, a peptide level output, and a protein level output,
#' whose names would be what is specified plus output levels.
#' For example, if the input string of \code{saveas} is "apir_output.xlsx", the saved PSM output would be "apir_output_psm.xlsx", the saved peptide output would be "apir_output_pepseq.xlsx", and the protein output would be "apir_output_pro.xlsx".
#' @param FDR the target FDR threshold
#' @param target_ls a named list of target PSM level output from search algorithms that need to combined.
#' Each element of the list is a data frame object.
#' The names should specify the search algorithms (see examples below).
#' We recommend users to input the complete target PSM level output by setting the FDR of search algorithms to be 1.
#' @param decoy_ls a list of decoy PSM level output from search algorithms that need to combined.
#' Each element of the list is a data frame object.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' We recommend users to input the complete decoy PSM level output by setting the FDR of search algorithms to be 1.
#' @param ifadjust a logical vector specifying if in the first round of identification q-value/pep thresholding should be applied, instead of APIR-adjust, to each search algorithm.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' We suggest users input a vector of \code{TRUE} if users are not sure the FDR control of individual search algorithms.
#' @param scoreColTitle a string vector specifying the column names of scores in \code{target_ls}.
#' Note that within a search algorithm, the column names of scores in the target output should match that in the decoy output.
#' We recommend users to input the column names of q-values or posterior error probabilities.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param fileColTitle a string vector specifying the column names of mass spectrum files in \code{target_ls}.
#' Note that within a search algorithm, the column names of mass spectrum files in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param scannumColTitle a string vector specifying the column names of scan numbers in \code{target_ls}.
#' Note that within a search algorithm, the column names of scan numbers in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param sequenceColTitle a string vector specifying the column names of peptide sequences in \code{target_ls}.
#' Note that within a search algorithm, the column names of peptide sequences in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param ifRecommendMasterProtein a logical specifying if a master protein should be aggregated using majority note.
#' If specified as \code{TRUE}, \code{masterproteinColTitle} also needs to be specified.
#' @param masterproteinColTitle a string vector specifying the column names of master proteins in \code{target_ls}.
#' Note that within a search algorithm, the column names of master proteins in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param ifRecommendModification a logical specifying if modifications should be aggregated.
#' If specified as \code{TRUE}, \code{modificationColTitle} also needs to be specified.
#' @param modificationColTitle a string vector specifying the column names of modifications in \code{target_ls}.
#' Note that within a search algorithm, the column names of modifications in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param staticModification a vector of static modification sites.
#' Could be one or more from "\code{C}"(Carbamidomethyl), "\code{K}"(TMT6plex), and "\code{N-term}"(TMT6plex).
#' The default is \code{c('C','K','N-term')}.
#' @param proteinPositionColTitle a string vector specifying the column names of modification positions in proteins in \code{target_ls}.
#' Note that within a search algorithm, the column names of master proteins in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param phospho_dataset a string specifying the location of the Phosphorylation_site_dataset downloaded from PhosphoSitePlus.
#' If not supplied, modifications sites are aggregated only based in frequency.
#' If supplied, modification sites are aggregated twice, once by frequency, once by both frequency and literature. See details.
#' @param organism a character string specifying the organism of your sample.
#' Could be any type available in the \code{ORGANISM} of \code{phospho_dataset}.
#' Only needed if phospho_dataset is not NULL.
#' @param ifAggregateAbundance a logical specifying if abundances should be aggregated.
#' If specified as \code{TRUE}, \code{abundanceColTitle} also needs to be specified.
#' @param abundanceColTitle a list of string vectors specifying the abundance column names in each search algorithm.
#' #' Note that within a search algorithm, the column names of abundances in the target output should match that in the decoy output.
#' Its length should be of the same as \code{target_ls}, and the order of search algorithms stored should match that in \code{target_ls}.
#' @param ncores the number of cores for parallel computing. The default is \code{parallel::detectCores()-1}.
#'
#'
#' @details
#' \code{apir} combines PSM level output from multiple search algorithms.
#' The following arguments, if specified, should have the same lengths as the number of search algorithms to be combined and should have the same order of search algorithms:
#' \code{target_ls}, \code{decoy_ls}, \code{ifadjust}, \code{scoreColTitle}, \code{scannumColTitle},
#' \code{sequenceColTitle}, \code{fileColTitle},\code{modificationColTitle}, \code{masterproteinColTitle},
#' \code{proteinPositionColTitle}, \code{abundanceColTitle}.
#'
#'
#'
#'
#'
#' @return \code{apir} outputs three excel files: one for the PSM level output, one for the peptide level, and another one for the protein level.
#' In its PSM level output, \code{apir} will append all columns that are in the target output from search algorithms.
#' @export
#' @importFrom parallel mclapply detectCores
#' @importFrom openxlsx write.xlsx read.xlsx
#' @import readr
#' @references
#' @author Yiling Chen, \email{yiling0210@ucla.edu}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
#'
#' @examples
#'
#' ### specifying arguments
#' target_ls = list(maxquant_target, msgf_target)
#' names(target_ls) = c('maxquant', 'msgf')
#' decoy_ls = list(maxquant_decoy, msgf_decoy)
#'
#' ifadjust = c(T, T)
#' scoreColTitle = c('PEP','QValue')
#' fileColTitle = c("Raw.file", 'SpectrumFile')
#' scannumColTitle = c('MS/MS.scan.number', "ScanNum")
#' sequenceColTitle = c('Sequence', 'Sequence')
#' masterproteinColTitle = c('Leading.razor.protein.(format.fixed)','Protein.(format.fixed)')
#' modificationColTitle = c("Modifications.ProteomeDiscoverer.Format","Modifications.ProteomeDiscoverer.Format")
#' proteinPositionColTitle = c("Positions.in.Proteins","Positions.in.Proteins")
#' abundanceColTitle = list(maxquant = c("Reporter.intensity.corrected.1", "Reporter.intensity.corrected.2","Reporter.intensity.corrected.3", "Reporter.intensity.corrected.4"),msgf = c('TMT126','TMT127','TMT128','TMT129'))
#'
#' ### no phosphosite based recommendation
#' apir(saveas = 'apir_output_wphospoplus.xlsx',
#' FDR = 0.05,
#' target_ls = target_ls,
#' decoy_ls = decoy_ls,
#' ifadjust = ifadjust,
#' scoreColTitle = scoreColTitle,
#' fileColTitle = fileColTitle,
#' scannumColTitle = scannumColTitle,
#' sequenceColTitle = sequenceColTitle,
#' ifRecommendMasterProtein = T,
#' masterproteinColTitle = masterproteinColTitle,
#' ifRecommendModification = T,
#' modificationColTitle = modificationColTitle,
#' staticModification = c('C','K','N-term'),
#' proteinPositionColTitle = proteinPositionColTitle,
#' organism = NULL,
#' phospho_dataset = NULL,
#' ifAggregateAbundance = T,
#' abundanceColTitle = abundanceColTitle)
#'
#' ###
#' organism = 'human'
#'
#' apir(saveas = 'apir_output.xlsx',
#' FDR = 0.05,
#' target_ls = target_ls,
#' decoy_ls = decoy_ls,
#' ifadjust = ifadjust,
#' scoreColTitle = scoreColTitle,
#' fileColTitle = fileColTitle,
#' scannumColTitle = scannumColTitle,
#' sequenceColTitle = sequenceColTitle,
#' ifRecommendMasterProtein = T,
#' masterproteinColTitle = masterproteinColTitle,
#' ifRecommendModification = T,
#' modificationColTitle = modificationColTitle,
#' staticModification = c('C','K','N-term'),
#' proteinPositionColTitle = proteinPositionColTitle,
#' organism = organism,
#' phospho_dataset = PhosphoSitePlus,
#' ifAggregateAbundance = T,
#' abundanceColTitle = abundanceColTitle)
#'
#'



apir = function(saveas,
                FDR,
                target_ls,
                decoy_ls,
                ifadjust,
                scoreColTitle,
                fileColTitle,
                scannumColTitle,
                sequenceColTitle,
                ifRecommendMasterProtein,
                masterproteinColTitle,
                ifRecommendModification,
                modificationColTitle,
                staticModification,
                proteinPositionColTitle,
                organism,
                phospho_dataset,
                ifAggregateAbundance,
                abundanceColTitle,
                ncores = parallel::detectCores() - 1){

  method_name = names(target_ls)
  n_method = length(target_ls)
  if(length(decoy_ls)!=n_method){
    stop('target and decoy data do not have the length')
  }
  if(length(ifadjust)!=n_method){
    stop('ifadjust does not have the same length as methods')
  }
  if(is.null(scoreColTitle) | is.null(scannumColTitle) | is.null(sequenceColTitle)){
    stop('variable names for q-values, scan number, and peptide sequences all need to be sepcified')
  }
  if(!length(scoreColTitle)%in% c(n_method,1) | !length(scannumColTitle)%in% c(n_method,1) | !length(sequenceColTitle)%in% c(n_method,1)){
    stop('Length of scoreColTitle, scannumColTitle, and sequenceColTitle have to be 1 or the same as target_ls')
  }
  if(ifAggregateAbundance & is.null(abundanceColTitle)){
    stop('When ifAggregateAbundance is TRUE, abundanceColTitle has to be specified')
  }
  if(ifRecommendModification & is.null(modificationColTitle)){
    stop('When ifRecommendModification is TRUE, modificationColTitle has to be specified')
  }
  if(!is.null(abundanceColTitle)){
    if(length(unique(sapply(abundanceColTitle,length)))>1){
      stop('Number of abundance variable is not consistent between methods!')
    }
  }
  if(!ifRecommendMasterProtein & ifRecommendModification){
    stop('ifRecommendMasterProtein has to be TRUE When ifRecommendModification is TRUE')
  }

  # object_to_max = match.arg(object_to_max, choices = c('psm','peptide','protein'), several.ok = F)
  scoreColTitle = rep(scoreColTitle, n_method/length(scoreColTitle))
  scannumColTitle = rep(scannumColTitle, n_method/length(scannumColTitle))
  sequenceColTitle = rep(sequenceColTitle, n_method/length(sequenceColTitle))


  target_ls_tot = target_ls
  target_ls_tot = lapply(1:n_method, function(i){
    x = target_ls_tot[[i]]
    query = paste( x[,fileColTitle[i]], x[,scannumColTitle[i]],sep = ':')
    match_id = paste( x[,fileColTitle[i]], x[,scannumColTitle[i]], x[,sequenceColTitle[i]],sep = ':')
    x$query = query
    x$match_id = match_id
    x
  })
  target_ls = lapply(1:n_method, function(i){
    x = target_ls_tot[[i]]
    ### modify scores such that large value means good, small value means bad
    scores = x[, scoreColTitle[i]]
    scores[scores ==0] = min(scores[scores>0])
    scores = -log10(scores)
    x = x[,c('query','match_id',sequenceColTitle[i])]
    x = cbind.data.frame(x, scores = scores, stringsAsFactors =F)
    x[!duplicated(x$match_id),]
  })
  decoy_ls = lapply(1:n_method, function(i){
    x = decoy_ls[[i]]
    query = paste( x[,fileColTitle[i]], x[,scannumColTitle[i]],sep = ':')
    match_id = paste( x[,fileColTitle[i]], x[,scannumColTitle[i]], x[,sequenceColTitle[i]],sep = ':')
    x$query = query
    x$match_id = match_id
    scores = x[, scoreColTitle[i]]
    scores[scores ==0] = min(scores[scores>0])
    scores = -log10(scores)
    x = x[,c('query','match_id',sequenceColTitle[i])]
    x = cbind.data.frame(x, scores = scores, stringsAsFactors =F)
    x[!duplicated(x$match_id),]  ### get rid of duplicated matches in decoy

  })
  ######## delete false matches from the decoy results #######
  decoy_ls = lapply(1:n_method, function(i){
    # print(i)
    dat_target = target_ls[[i]]
    dat_decoy = decoy_ls[[i]]
    match_id = dat_target$match_id
    match_decoy = dat_decoy$match_id
    false_decoy =  match_decoy %in% match_id
    # print(sum(false_decoy))
    dat_decoy[!false_decoy,]
  })

  ###############
  seq_pro_tot = lapply(1:n_method, function(i){
    x = target_ls_tot[[i]]
    temp = x[, c(sequenceColTitle[i],masterproteinColTitle[i])]
    names(temp) = c('seq','pro')
    temp
  })
  seq_pro_tot = Reduce('rbind.data.frame', seq_pro_tot)
  dim(seq_pro_tot)
  seq_pro_tot = split(seq_pro_tot$pro, f = factor(seq_pro_tot$seq))
  seq_pro_tot = lapply(seq_pro_tot, function(x){
    x_new = unique(unlist(strsplit(x, split='; ')))
    x_new[!is.na(x_new)]
    return(x_new)
  })


  method_added = integer(0)
  method_left = setdiff(1:n_method,  method_added)

  final_disc = list()
  match_examined = list()

  ########## find the first set ############
    psm_bym = lapply(method_left, function(i_method){
      # print(i_method)
      if( ifadjust[i_method] ){
        # print('adjust')
        dat_target = target_ls[[i_method]]
        dat_target = dat_target[! dat_target$match_id %in% unlist(match_examined),,drop = F ]
        dat_decoy = decoy_ls[[i_method]]
        if(nrow(dat_target) >0){
          psm_i = apir_select(FDR = FDR,
                              target_match = dat_target$match_id,
                              target_query = dat_target$query,
                              signal = dat_target$scores,
                              decoy_match = dat_decoy$match_id,
                              decoy_query = dat_decoy$query,
                              control = dat_decoy$scores,
                              ncores = ncores)
        }else{
          psm_i = character(0)
        }

        return(psm_i)

        }else{
          # print('not adjust')

          x = target_ls_tot[[i_method]]
          unique(x$match_id[x[,scoreColTitle[i_method]] <= FDR])
      }
      # print(i_method)

    })


  which2add = find_optimal_output( object_to_max = 'peptide',psm_bym)
  match_examined = append(match_examined, list(target_ls[[method_left[which2add]]]$match_id))
  method_added = c(method_added, method_left[which2add])
  method_left = setdiff(1:n_method,  method_added)
  final_disc = append(final_disc, list(psm_bym[[which2add]]))


  ################### add psms from other search engines
  while(length(method_left)>0){
    # print(paste0('have found ',length(unlist(final_disc)),' psm'))
    psm_bym = lapply(method_left, function(i_method){
      # print(i_method)
      dat_target = target_ls[[i_method]]
      dat_target = dat_target[! dat_target$match_id %in% unlist(match_examined),,drop = F ]
      dat_decoy = decoy_ls[[i_method]]
      if(nrow(dat_target) >0){
        psm_i = apir_select(FDR = FDR,
                            target_match = dat_target$match_id,
                            target_query = dat_target$query,
                            signal = dat_target$scores,
                            decoy_match = dat_decoy$match_id,
                            decoy_query = dat_decoy$query,
                            control = dat_decoy$scores,
                            ncores = ncores)
      }else{
        psm_i = character(0)
      }

      return(psm_i)
    })

      if( all(sapply(psm_bym,length) == 0 )){
        # print('The remaining search engines all found 0 psms')
        break
      }

      which2add = find_optimal_output( object_to_max = 'peptide',psm_bym)
      match_examined = append(match_examined, list(target_ls[[method_left[which2add]]]$match_id))
      # print(paste0('adding method ', method_left[which2add]))
      method_added = c(method_added, method_left[which2add])
      method_left = setdiff(1:n_method,  method_added)
      # print(paste0('method left to examine ', paste0(method_left, collapse = ' ')))
      final_disc = append(final_disc, list(psm_bym[[which2add]]))
    }



  # Reduce('intersect',final_disc)


  # method_added = method_name[method_added]
  # method_added = rep(method_added, sapply(final_disc, length))
  final_disc = unlist(final_disc)
  searchEngine = sapply(target_ls, function(x){
    (final_disc %in% x$match_id)
  })
  searchEngine = sapply(1:nrow(searchEngine), function(i){
    paste0(method_name[searchEngine[i,]],collapse = ';')
  })
  id_split = lapply(final_disc, strsplit, split = ':')
  rawfile = sapply(id_split, function(x){unlist(x)[1]})
  scannum = sapply(id_split, function(x){unlist(x)[2]})
  peptideseq = sapply(id_split, function(x){unlist(x)[3]})
  pro_tot = seq_pro_tot[match(peptideseq, names(seq_pro_tot))]
  pro_tot = sapply(pro_tot,function(x){
    paste0(x, collapse = ';')
  })

  remove(id_split, match_examined)

  ########################################

    print('Performing protein inference and aggregating abundances...')
    if(ifAggregateAbundance){
      abundance_methods = lapply(1:n_method, function(i){
        # print(i)
        x = target_ls_tot[[i]][,  abundanceColTitle[[i]]]
        x[match(final_disc, target_ls_tot[[i]]$match_id),]
      })

      abundance_aggregated = data.frame(colnormalize_abundance(abundance_methods))
      names(abundance_aggregated) = paste0('aggregatedAbundance.channel',1:length(abundanceColTitle[[1]]))
    }
    if(ifRecommendMasterProtein){
      masterproteins = lapply(1:n_method, function(i){
        # print(i)
        x = target_ls_tot[[i]][, masterproteinColTitle[[i]]]
        x[match(final_disc, target_ls_tot[[i]]$match_id)]
      })
      masterproteins = as.data.frame(masterproteins, stringsAsFactors = F)
      masterprotein_recommended = recommend_masterproteins(masterproteins)

      if(!is.null(proteinPositionColTitle)){
        seqprositionsinprotein = lapply(1:n_method, function(i){
          x = target_ls_tot[[i]][, proteinPositionColTitle[i]]
          x[match(final_disc, target_ls_tot[[i]]$match_id)]
        })
        seqprositionsinprotein = as.data.frame(seqprositionsinprotein)
        seqposition_recommended = recommend_seqposition(seqprositionsinprotein, masterprotein_recommended )

      }


    }

    if(ifRecommendModification){

      modifications_methods = lapply(1:n_method, function(i){
        # print(i)
        x = target_ls_tot[[i]][, modificationColTitle[i]]
        x[match(final_disc, target_ls_tot[[i]]$match_id)]
      })
      modifications_methods = as.data.frame(modifications_methods, stringsAsFactors = F)
      ############ if phospho_dataset is supplied, we will have recommend using that
      ############ otherwise, just majority vote
      modifications_recommended = recommend_modifications(method_name,
                                                        modifications_methods,
                                                        masterprotein_recommended,
                                                        phospho_dataset,
                                                        organism = organism,
                                                        staticModification = staticModification)
    }
    if(ifRecommendModification & ifRecommendMasterProtein){
      recommendmod_colname = ifelse(is.null(phospho_dataset),"modifications_recommended","modifications_recommended_phosphoSitePlus" )
      modifications_in_proteins = find_modifications_in_protein(seqposition_recommended = seqposition_recommended,
                                                                modifications_recommended = modifications_recommended[,recommendmod_colname])
    }



    print('Formatting the final output...')

    qval_methods = lapply(1:n_method, function(i){
      x = target_ls_tot[[i]][,scoreColTitle[i]]
      x[match(final_disc, target_ls_tot[[i]]$match_id)]
    })
    names(qval_methods) = paste0(method_name,sep = '.', 'q-values')

  if(!is.null(saveas)){

    modifications_methods = lapply(1:n_method, function(i){
      x = target_ls_tot[[i]][, modificationColTitle[i]]
      x[match(final_disc, target_ls_tot[[i]]$match_id)]
    })
    names(modifications_methods) = paste0(method_name,sep = '.', 'Modifications')

    out_core = cbind.data.frame(scannum = scannum,
                           rawfile = rawfile,
                           sequence = peptideseq,
                           searchEngine = searchEngine,
                           # qval_methods,
                           # masterprotein_union = pro_tot,
                           # modifications_methods,
                           stringsAsFactors = F
    )


    if(ifRecommendMasterProtein){
      out_core = cbind.data.frame(out_core,
                                  masterprotein_recommended = masterprotein_recommended,
                                  position_in_protein_recommended = seqposition_recommended, stringsAsFactors = F)
    }
    if(ifRecommendModification){
      out_core = cbind.data.frame(out_core, modifications_recommended, stringsAsFactors = F)
    }
    if(ifRecommendModification & ifRecommendMasterProtein){
      out_core = cbind.data.frame(out_core, modifications_in_proteins = modifications_in_proteins)
    }
    if(ifAggregateAbundance){
      out_core = cbind.data.frame(out_core, abundance_aggregated, stringsAsFactors = F)
    }


    target_ls_unused = lapply(1:n_method, function(i){
      x = target_ls_tot[[i]]
      x = x[match(final_disc, target_ls_tot[[i]]$match_id),
            !colnames(x) %in% c(scannumColTitle[i],fileColTitle[i],sequenceColTitle[i] )]
      names(x) = paste0(method_name[i], sep = '.',colnames(x))
      return(x)
    })
    out_psm = cbind.data.frame(out_core, as.data.frame(target_ls_unused, stringsAsFactors = F), stringsAsFactors = F)
    saveas_psm = gsub(x= saveas,pattern = '.xlsx', replacement = '_psm.xlsx' )
    write.xlsx(out_psm, file = saveas_psm)

    ######## start from here to add protein modifications

    ################## convert to peptide level output ##################
    out_pepseq = out_core
    pepseqmod = paste0(out_pepseq$sequence,'+',out_pepseq[,  recommendmod_colname ])

    out_pepseq$spectrumID = paste0(out_pepseq$rawfile, ":",out_pepseq$scannum)
    out_pepseq$scannum = NULL
    out_pepseq$rawfile = NULL
    out_pepseq = split.data.frame(out_pepseq, f = pepseqmod)
    colname_outpepseq = colnames(out_pepseq[[1]])
    out_pepseq = lapply(1:length(out_pepseq), function(i_x){
      # print(i_x)
      x_i = out_pepseq[[i_x]]

      modifications_in_proteins =  paste0(unique(unlist(strsplit(x_i[, "modifications_in_proteins"], split = '; '))), collapse = '; ')

      sequence = unique(x_i$sequence)

      spectrumID = paste0(unique(x_i[, 'spectrumID']), collapse = ';')

      searchEngine = paste0(unique(unlist(strsplit(x_i$searchEngine, split =';'))), collapse = ';')

      modifications_recommended = unique(x_i[, recommendmod_colname])

      abundance_aggregated = matrix(colSums(x_i[,grep(pattern = 'ggregatedAbundance.',x = colname_outpepseq)],na.rm = T),nrow =1)
      abundance_aggregated = as.data.frame(abundance_aggregated)
      abundance_aggregated = t(apply(abundance_aggregated, 1, function(v){
        # mean(unlist(
        v/sum(v,na.rm = T)*100*ncol(abundance_aggregated)
        # ))
      }))
      colnames(abundance_aggregated) =  colname_outpepseq[grep(pattern = 'ggregatedAbundance.',x = colname_outpepseq)]

      masterprotein_recommended = unique(x_i[,'masterprotein_recommended'])
      masterprotein_recommended = masterprotein_recommended[!is.na(masterprotein_recommended)]
      masterprotein_recommended = paste0(masterprotein_recommended,collapse =";")

      # modifications_recommended_phosphoSitePlus = unique(x_i[, 'modifications_recommended_phosphoSitePlus'])
      # modifications_recommended_phosphoSitePlus = paste0(modifications_recommended_phosphoSitePlus, collapse = ';')


      positions_rec = unique(x_i[, "position_in_protein_recommended"])
      positions_rec = paste0(positions_rec, collapse =";")


      x_i_new = cbind.data.frame(sequence = sequence,
                       modifications_recommended = modifications_recommended,
                       spectrumID = spectrumID,
                       searchEngine = searchEngine,
                       masterprotein_recommended,
                       position_in_protein_recommended = positions_rec,
                       modifications_in_proteins = modifications_in_proteins,
                       abundance_aggregated ,
                      stringsAsFactors = F
                       )
      return(x_i_new)
    })
    out_pepseq = Reduce('rbind.data.frame', out_pepseq)
    saveas_pepseq = gsub(x= saveas,pattern = '.xlsx', replacement = '_pepseq.xlsx' )
    write.xlsx(out_pepseq, file = saveas_pepseq)


    ################## convert to protein level output ##################
    out_pro = out_core
    masterpro = out_pro$masterprotein_recommended
    out_pro$spectrumID = paste0(out_pro$rawfile, ":",out_pro$scannum)
    out_pro$scannum = NULL
    out_pro$rawfile = NULL
    out_pro$sequence = NULL
    out_pro$position_in_protein_recommended = NULL
    out_pro$modifications_recommended = NULL
    out_pro$modifications_recommended_phosphoSitePlus = NULL
    # out_pro$position_in_protein_recommended = NULL
    out_pro = split.data.frame(out_pro, f = masterpro)
    colname_outpro = colnames(out_pro[[1]])

    out_pro = lapply(1:length(out_pro), function(i_x){
      # print(i_x)
      x_i = out_pro[[i_x]]
      masterprotein_recommended = x_i$masterprotein_recommended
      masterprotein_recommended = unique(masterprotein_recommended)

      modifications_in_proteins = x_i$modifications_in_proteins
      modifications_in_proteins = unique(modifications_in_proteins)
      modifications_in_proteins =  modifications_in_proteins[!is.na(modifications_in_proteins)]
      if(length(modifications_in_proteins) == 0){
        modifications_in_proteins_formatted  = NA
      }else{
        modifications_in_proteins_formatted = format_modifications_in_proteins(modifications_in_proteins)

      }

      spectrumID = x_i$spectrumID
      spectrumID = paste0(unique(spectrumID), collapse = '; ')

      searchEngine = paste0(unique(unlist(strsplit(x_i$searchEngine, split =';'))), collapse = ';')

      aggregated_abundance = x_i[, grepl('Abundance', colname_outpro)]
      aggregated_abundance = colSums(aggregated_abundance)
      aggregated_abundance = aggregated_abundance/sum( aggregated_abundance, na.rm =T)*100*length(aggregated_abundance)
      aggregated_abundance = data.frame(matrix(aggregated_abundance, nrow = 1))
      colnames(aggregated_abundance) = colname_outpro[grepl('Abundance', colname_outpro)]
      x_i_new = cbind.data.frame(masterprotein_recommended = masterprotein_recommended,
                                 modifications_in_proteins = modifications_in_proteins_formatted,
                                 spectrumID = spectrumID,
                                 searchEngine = searchEngine,
                                 aggregated_abundance, stringsAsFactors = F)
      return(x_i_new)
    })
    out_pro = Reduce('rbind.data.frame', out_pro)
    saveas_pro = gsub(x= saveas,pattern = '.xlsx', replacement = '_pro.xlsx' )
    write.xlsx(out_pro, file = saveas_pro)
  }

    # re = list(psm = final_disc,
    #           pepseq = peptideseq,
    #           masterpro_tot = as.vector(pro_tot),
    #           masterpro_recommend = masterprotein_recommended)
    # if(ifRecommendModification){
    #   re = append(re, list(modifications_recommended_phosphoSitePlus = modifications_recommended[,2],
    #                        modifications_recommended = modifications_recommended[,1]))
    # }
    #
    # return(re)
    return(NULL)
}
