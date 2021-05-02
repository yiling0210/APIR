library(parallel)
library(openxlsx)
library(readr)
# parse_number = function (x, na = c("", "NA"), locale = default_locale(), trim_ws = TRUE){
#   parse_vector(x, col_number(), na = na, locale = locale, trim_ws = trim_ws)
# }

clipper = function(signal, control, FDR,  ncores = detectCores() - 1 ){


  contrastScore = pmax(signal, control)* sign(signal -  control)

  re = clipper_w_c(contrastScore = contrastScore, FDR= FDR, ncores = ncores)

  return(re)
}
# mod = c("S676(Phospho);S686(Phospho);M690(Oxidation)")
format_modifications_in_proteins = function(modifications_in_proteins){
  prot = sapply(strsplit(modifications_in_proteins, split = ' '),'[[',2)
  mod = sapply(strsplit(modifications_in_proteins, split = ' '),'[[',2)
  mod = unlist(strsplit(mod,';'))
  mod = unlist(strsplit(mod, '/'))
  mod = unique(mod)
  mod = mod[grepl(pattern = 'Phospho',mod) | grepl(pattern = 'Oxidation',mod)]
  mod = mod[order(parse_number(mod), decreasing = F)]
  if(length(mod) == 0){
    re = NA
  }else{
    mod = unlist(strsplit(mod, split = '[()]'))
    mod_pos = mod[seq(1, length(mod)-1, by = 2)]
    mod_type = mod[-seq(1, length(mod)-1, by = 2)]
    re = split(mod_pos, as.factor(mod_type))
    counts = sapply(re, length)
    re = lapply(re, paste0, collapse = '; ')
    re = paste0(counts, '\u00D7',names(re),' [', unlist(re),']' )
    re = paste0(re, collapse = '; ')
  }

  return(re)
}

clipper_w_c = function(contrastScore, FDR, ncores = detectCores() -1 ){
  contrastScore[is.na(contrastScore)] = 0
  # print(hist(contrastScore))
  c_abs = abs(contrastScore[contrastScore != 0])
  c_abs  = sort(unique(c_abs))

  emp_fdp = mclapply(c_abs, function(t){
    (1 + sum(contrastScore <= -t))/ sum(contrastScore >= t)
  },mc.cores = ncores)
  emp_fdp = unlist(emp_fdp)

  idx = which(emp_fdp <= FDR)

  if(length(idx) == 0){
    re = list(
      discovery = which(1<0),
      # contrastScore = contrastScore,
      thre =  which(1<0))
  }else{
    thre = c_abs[min(idx)]
    re = list(
      discovery = which( contrastScore>= thre),
      thre = thre
    )

  }


  return(re)
}


pooled_one_sample = function(signal, control,FDR, ncores){
  pval = mclapply(signal, function(s){
    mean(control > s)
  },mc.cores = ncores)
  pval = unlist(pval)
  pval_adj = p.adjust(pval, method = 'BH')
  discovery = which(pval_adj <= FDR)
  return(discovery)
}

# FDR = FDR
#
# target_match = dat_target$match_id[idx2select]
# target_query = dat_target$query[idx2select]
# signal = dat_target$scores[idx2select]
# decoy_match = dat_decoy$match_id
# decoy_query = dat_decoy$query
# control = dat_decoy$scores
# ncores = ncores
# pair_thre = 0.4
# FDR = FDR
# target_match = dat_target$match_id
# target_query = dat_target$query
# signal = dat_target$scores
# decoy_match = dat_decoy$match_id
# decoy_query = dat_decoy$query
# control = dat_decoy$scores
# ncores = ncores
apir_select = function(FDR, target_match, target_query, signal, decoy_match, decoy_query, control, pair_thre = 0.4, ncores ){
  perct_covered = mean(target_query %in% decoy_query)
  perct_covered
  if(perct_covered > pair_thre){
    #### use clipper
    shared_query = target_query[target_query %in% decoy_query]
    shared_match = target_match[target_query %in% decoy_query]
    x = signal[match(shared_match, target_match)]
    y = control[match(shared_query, decoy_query)]
    re_clipper_m = clipper(signal = x, control = y,
                           FDR = FDR, ncores = ncores )
    matches_clipper_m = shared_match[re_clipper_m$discovery]
    return(unique(matches_clipper_m))

  }else{
    #### use pooled approach
    re_pooled_m = pooled_one_sample(signal = signal, control = control,
                                    FDR = FDR, ncores = ncores )
    return(unique(target_match[re_pooled_m]))
  }


}


find_optimal_output = function(object_to_max, matches_bym){
  if(object_to_max == 'psm'){
    base_m = which.max(sapply(matches_bym, length))
  }
  if(object_to_max %in% c('peptide')){
    id_split_bym = lapply(matches_bym, strsplit, split = ':')

    peptideseq_bym = sapply(id_split_bym, function(x){
      unique(sapply(x,'[',3))
    })

    if(object_to_max  == 'peptide'){
      base_m = which.max(sapply(peptideseq_bym,length))
    }


  }
  return(base_m)
}

colnormalize_abundance = function(abundance_methods){

  n_r = nrow(abundance_methods[[1]])
  n_c = ncol(abundance_methods[[1]])
  abundance_methods = sapply(abundance_methods, function(x){
    as.matrix(x, ncol = 1, byrow = F)
  })
  ###############
  avg_abundance = rowMeans(abundance_methods, na.rm = T)
  avg_abundance[is.nan(avg_abundance)] = NA
  avg_abundance = matrix(avg_abundance, ncol = n_c, nrow = n_r, byrow = F )

  f = (10^6)/apply(avg_abundance,2,sum,na.rm = T)
  norm_abundance = t(apply(avg_abundance, 1, "*",f))
  return(norm_abundance)
  # abundance_scaled = apply(norm_abundance , 1, function(x){
  #   x/mean(x)*100
  # })
  # abundance_scaled[is.nan(abundance_scaled)] = NA
  # return(t(abundance_scaled))
}

recommend_masterproteins = function(masterproteins, ncores ){
  temp = mclapply(1:nrow(masterproteins), function(i){
    x = masterproteins[i,]
    x = as.character(x)
    if(all(is.na(x))){
      return(NA)
    }else{
      x = unlist(strsplit(x, '; '))
      tb = table(x)
      return(names(tb)[which.max(tb)])
    }
  }, mc.cores = ncores)

  temp = unlist(temp)

  return(temp)
}

recommend_seqposition = function(seqprositionsinprotein, masterprotein_recommended ){
  seq_position = mapply(function(masterpro, i){
    # print(i)
    positions = unlist(seqprositionsinprotein[i,])
    if(all(is.na(positions))){
      return(NA)
    }else{
      cand = positions[grep(masterpro, (positions))]
      cand = strsplit(cand, '; ')
      cand = sapply(cand, function(cand_i){
        cand_i = cand_i[grep(masterpro, cand_i)]
        paste0(cand_i, collapse = '; ')
      })
      tb = table(cand)
      paste0(names(tb)[which(tb == max(tb))], collapse = ';')
    }

  }, masterpro = masterprotein_recommended, i = 1:nrow(seqprositionsinprotein), SIMPLIFY = T)
  return(seq_position)
}

find_most_cited = function(phospho, idx, cand){
  citations = phospho[, c('MS_LIT','MS_CST')]
  citations[is.na(citations)] = 0
  # citations$`CST_CAT.` = as.numeric(citations$`CST_CAT.`)
  citations = citations[idx[!is.na(idx)],]
  rk = NA
  j = 1
  while(is.na(rk) & j < 3){
    if(length(unique(citations[,j])) == 1){
      j = j + 1
    }else{
      rk = which.max(citations[,j])
    }

  }
  if(is.na(rk)){
    return(paste0(sort(cand[!is.na(idx)]),collapse = '/'))

  }else{
    return(cand[!is.na(idx)][rk])
  }



}

collapse2modifications = function(mod_ls){
  type_mod = substr(mod_ls, start = 1, stop = 1)
  if( all(c('C','S') %in% type_mod) ){
    out = mod_ls[which(type_mod == 'S')]
  }else if( all(c('N','S') %in% type_mod) ){
    out = mod_ls[which(type_mod == 'S')]
  }else if( all(c('K','S') %in% type_mod) ){
    out = mod_ls[which(type_mod == 'S')]
  }else{
    out = paste0(mod_ls, collapse = '/')
  }
  return(out)
}

recommend_modifications = function(method_name,
                                   modifications_methods,
                                   masterprotein_recommended,
                                   phospho_dataset,
                                   organism,
                                   staticModification,
                                   ncores
){

  staticModification = substr(staticModification, 1,1)
  #### find the recommended masterprotein

  n_matches = nrow(modifications_methods)
  n_method = ncol(modifications_methods)
  modifications_methods = sapply(1:n_method, function(i){
    gsub(modifications_methods[,i], replacement = '', pattern = '\\(Prot\\)')
  })


  recommended_modifications = mclapply(1:n_matches,function(i){

    x = as.character(modifications_methods[i,])
    h = strsplit(x, split ='; ')
    ##### extract static modifications #####
    staticmod = h[!is.na(h)][[1]]
    staticmod = staticmod[substr(staticmod, 1,1) %in% staticModification]

    ##### extract varying modifications #####
    varyingmod = lapply(h, function(h_i){
      h_i[!h_i %in% staticmod]
    })
    n_modifications = max(sapply(varyingmod, length))

    tb = table(unlist(varyingmod))

    varyingmod = lapply(varyingmod, function(temp){
      temp[order(tb[match( temp, names(tb))], decreasing = T)]
    }) ## order varying modifications by frequency and category

    ### recommended by majority rule only
    i_position = 1
    re_nositeplus = rep(NA,0)
    while(i_position <= n_modifications){

      cand = sapply(varyingmod, '[',i_position)

      cand_simplified = sapply(cand, strsplit, split = '[()]')
      cand_simplified = sapply(cand_simplified, '[',1)
      tb = table(cand)
      if(sum(tb == max(tb)) == 1){
        re_nositeplus = c(re_nositeplus, names(which.max(tb)))
      }
      if(sum(tb == max(tb)) > 1){
        mod_ls = sort(names(tb)[tb == max(tb)])
        re_nositeplus = c(re_nositeplus, collapse2modifications(mod_ls))

      }
      i_position = i_position + 1
    }
    re_nositeplus = c(staticmod, re_nositeplus)
    positions = parse_number(re_nositeplus[-1])
    re_nositeplus = c(re_nositeplus[1],re_nositeplus[-1][order(positions)] )
    re_nositeplus = paste(re_nositeplus, collapse = ";")


    if(!is.null(phospho_dataset)){
      i_position = 1
      re_siteplus = rep(NA,0)
      while(i_position <= n_modifications){
        # print(paste0(i,'-th row; ', i_position, '-th position'))

        cand = sapply(varyingmod, '[',i_position)

        cand_simplified = sapply(cand, strsplit, split = '[()]')
        cand_simplified = sapply(cand_simplified, '[',1)
        tb = table(cand)
        if(sum(tb == max(tb)) == 1){
          re_siteplus = c(re_siteplus, names(which.max(tb)))
        }
        if(sum(tb == max(tb)) > 1){
          masterprotein_i = masterprotein_recommended[i]
          phospho =  phospho_dataset[phospho_dataset$ACC_ID == masterprotein_i & phospho_dataset$ORGANISM == organism,]
          modifications_masterprotein_i = phospho$MOD_RSD
          modifications_masterprotein_i = strsplit(modifications_masterprotein_i, split = "-|\\s")
          modifications_masterprotein_i = sapply(modifications_masterprotein_i, '[',1)

          idx = match(names(tb)[tb == max(tb)],modifications_masterprotein_i)
          # idx = idx[!is.na(idx)]
          if(all(is.na(idx))){
            re_siteplus = c(re_siteplus, paste0(names(tb)[tb == max(tb)],collapse  = '/'))
          }
          if(sum(!is.na(idx)) == 1){
            re_siteplus = c(re_siteplus, unique(cand[!is.na(idx)]))
          }
          if(sum(!is.na(idx)) > 1){
            # find_most_cited(phospho, idx, cand)
            re_siteplus = c(re_siteplus,find_most_cited(phospho, idx, cand))
          }

        }
        i_position = i_position + 1
      }

      re_siteplus = c(staticmod, re_siteplus)
      positions = parse_number(re_siteplus[-1])
      re_siteplus = c(re_siteplus[1],re_siteplus[-1][order(positions)] )
      re_siteplus = paste(re_siteplus, collapse = ";")
    }


    if(!is.null(phospho_dataset)){
      return(c(re_nositeplus, re_siteplus))
    }else{
      return(re_nositeplus)
    }


  }, mc.cores = ncores)

  recommended_modifications = matrix(unlist(recommended_modifications), byrow = T, nrow = n_matches)

  # if(is.null(dim(recommended_modifications))){
  #   recommended_modifications = matrix(recommended_modifications, ncol = 1)
  # }else{
  #   recommended_modifications = t(recommended_modifications)
  # }

  colnames(recommended_modifications) = if(is.null(phospho_dataset)){
    'modifications_recommended'
  }else{
    c('modifications_recommended','modifications_recommended_PSP')
  }
  return( recommended_modifications)

}

# mod = m
# pos = p
find_modifications_in_protein_byPSM = function(mod, pos){
  prot = strsplit(pos, split = ' ')[[1]][1]  # extract protein from position
  pos = strsplit(pos, split = ' ')[[1]][2] # extract peptide position from position
  st_pos = parse_number(pos)
  mod_type = regmatches(mod, gregexpr("(?=\\().*?(?<=\\))", mod, perl=T))[[1]] # extract mod type with parenthesis

  mod = unlist(strsplit(mod, split = '[()]')) # separate mod type from sites
  mod_pos = mod[seq(1, length(mod)-1, by = 2)] # extract sites by extracting the odd number position
  i = 1
  mod_pos_new = rep(NA, length(mod_pos))
  while(i <= length(mod_pos)){
    mod_pos_i = mod_pos[i]
    if(!grepl(pattern = "[[:digit:]]+", mod_pos_i)){
      mod_pos_new[i] = mod_pos_i
      i = i + 1
    }else{
      num = unlist(regmatches(mod_pos_i, gregexpr("[[:digit:]]+", mod_pos_i)))
      mod_pos_new[i] = gsub(pattern = num, replacement = as.character(as.numeric(num) + st_pos - 1), mod_pos_i)
      i = i + 1
    }
  }
  mod_ls = c(mod_pos_new, mod_type)
  mod_ls = mod_ls[rep(1:(length(mod)/2), each = 2) + c(0,length(mod)/2 )]
  re = paste0(mod_ls, collapse = '')
  re = paste0(prot, ' ', re)
  return(re)
}

# seqposition_recommended = seqposition_recommended
# modifications_recommended = modifications_recommended[,recommendmod_colname]
find_modifications_in_protein = function(seqposition_recommended, modifications_recommended, ncores){
  n_match = length(seqposition_recommended)
  modifications_in_proteins = mclapply(1:n_match, function(i){

    m = modifications_recommended[i]
    m = paste0(unlist(strsplit(m, split = ';'))[-1], collapse = ';')

    p = seqposition_recommended[i]
    if(m == ""){
      return(NA)
    }else if(!grepl(';', p)){
      re = find_modifications_in_protein_byPSM(mod = m, pos = p)
    }else{
      p_ls = unlist(strsplit(p, split = ';'))
      p_ls_new = sapply(p_ls, function(p_i){
        find_modifications_in_protein_byPSM(mod= m, pos = p_i)
      })
      re = paste0(p_ls_new, collapse = '; ')
    }

    return(re)
  }, mc.cores = ncores)
  modifications_in_proteins = unlist(modifications_in_proteins)
  return(modifications_in_proteins)
}
