def read_tar_inner_filenames(tar_path):
    """
    读取 TAR 压缩包内所有文件/文件夹的名称
    :param tar_path: TAR 文件路径（支持 .tar/.tar.gz/.tgz/.tar.bz2 等）
    :return: 列表，包含压缩包内所有文件/文件夹名称；失败返回空列表
    """
    
    if not os.path.exists(tar_path):
        print(f"error：file {tar_path} didn't exist！")
        return []
    
    try:
        with tarfile.open(tar_path, mode='r') as tf:
            
            inner_filenames = tf.getnames()
            
            print(f"TAR 包 {tar_path} 内的文件/文件夹列表：")
            for idx, name in enumerate(inner_filenames, 1):
                print(f"{idx}. {name}")
        
        return inner_filenames
    
    except tarfile.TarError as e:
        print(f"错误：TAR 文件格式异常/损坏 → {str(e)}")
        return []
    except PermissionError:
        print(f"错误：无权限读取 {tar_path}")
        return []
    except Exception as e:
        print(f"读取失败：{str(e)}")
        return []
    
    adata.obs_names = [f"{sample_name}_{cell_id}" for cell_id in adata.obs_names]
    
    if not adata.var_names.is_unique:
        dup_genes = adata.var_names[adata.var_names.duplicated()].unique()
        new_var_names = []
        gene_count = {}
        for gene in adata.var_names:
            if gene in gene_count:
                gene_count[gene] += 1
                new_var_names.append(f"{gene}_{gene_count[gene]}")
            else:
                gene_count[gene] = 0
                new_var_names.append(gene)
        adata.var_names = new_var_names
        
    return adata

stages = ['normal', 'N', 'HSIL', 'SCC', 'ADC']

if not os.listdir(extract_dir):
    with tarfile.open(tar_path, 'r') as tar:
        tar.extractall(path=extract_dir)

all_processed_data = {}

for group_name, gsm_dict in sample_mapping.items():
    
    print(f"\n{'~' * 10} processing {group_name} {'~' * 10} ")
    
    group_adata_list = []
    
    for gsm_prefix, sample_name in gsm_dict.items():
        matrix_path = os.path.join(extract_dir, f"{gsm_prefix}_{sample_name}.matrix.mtx.gz")
        barcodes_path = os.path.join(extract_dir, f"{gsm_prefix}_{sample_name}.barcodes.tsv.gz")
        features_path = os.path.join(extract_dir, f"{gsm_prefix}_{sample_name}.features.tsv.gz")
        
        if not all([os.path.exists(f) for f in [matrix_path, barcodes_path, features_path]]):
            print(f"{sample_name}：文件不完整")
            exit()
        
        try:
            with gzip.open(matrix_path, 'rb') as f:
                expr_matrix = mmread(f).tocsr()
            
            with gzip.open(barcodes_path, 'rt') as f:
                barcodes = [line.strip() for line in f.readlines()]
            
            with gzip.open(features_path, 'rt') as f:
                features = []
                for line in f.readlines():
                    line = line.strip()
                    if line:
                        features.append(line.split('\t')[1])
            
            adata = ad.AnnData(X=expr_matrix.T)
            adata.obs_names = barcodes[:adata.n_obs]
            adata.var_names = features[:adata.n_vars]
            
            adata = fix_duplicate_labels(adata, sample_name)
            
            adata.obs['sample'] = sample_name
            adata.obs['group'] = group_name
            
            print(f" {sample_name}：{adata.n_obs}细胞 × {adata.n_vars}基因")
            group_adata_list.append(adata)
        
        except Exception as e:
            
            print(f"处理 {sample_name} 失败：{str(e)}")
            continue
    
    if group_adata_list:
        group_adata = ad.concat(
            group_adata_list, 
            label='sample', 
            join='outer',
            fill_value=0,
            index_unique='-'
        )
        
        sc.pp.calculate_qc_metrics(group_adata, percent_top=None, log1p=False, inplace=True)
        
        group_adata.var['mt'] = group_adata.var_names.str.contains('^MT-|^mt-', regex=True)
        # 更全面的核糖体基因匹配
        group_adata.var['ribo'] = group_adata.var_names.str.contains('^RPS|^RPL|^Rps|^Rpl', regex=True)
        
        sc.pp.calculate_qc_metrics(
            group_adata, 
            qc_vars=['mt', 'ribo'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        
        group_adata = group_adata[
            (group_adata.obs['n_genes_by_counts'] > QC_PARAMS['min_genes']) &
            (group_adata.obs['n_genes_by_counts'] < QC_PARAMS['max_genes']) &
            (group_adata.obs['pct_counts_mt'] < QC_PARAMS['max_mt_pct']), :
        ]
        
        sc.pp.filter_genes(group_adata, min_cells=QC_PARAMS['min_cells'])

        #去除基因
        group_adata = group_adata[:, ~(group_adata.var['mt'] | group_adata.var['ribo'])].copy()
        
        print(f"QC后 {group_name}：{group_adata.n_obs}细胞 × {group_adata.n_vars}基因")
        
        expr_df = pd.DataFrame(
            group_adata.X.toarray() if hasattr(group_adata.X, 'toarray') else group_adata.X,
            index=group_adata.obs_names,
            columns=group_adata.var_names
        ).T
        
        csv_path = os.path.join(output_dir, f"{group_name}_raw_counts.csv")
        h5ad_path = os.path.join(output_dir, f"{group_name}_processed.h5ad")
        expr_df.to_csv(csv_path, index=True)
        group_adata.write_h5ad(h5ad_path)
        
        all_processed_data[group_name] = expr_df
        print(f"保存完成：{csv_path}")
    else:
        print(f"{group_name} 组无有效数据！")

    adata = group_adata if group_adata_list else ad
    available_genes = set(adata.var_names)
    cd4_exists = 'CD4' in available_genes
    cd8a_exists = 'CD8A' in available_genes
    cd8b_exists = 'CD8B' in available_genes
    
    print(f"  基因存在性检查: CD4={cd4_exists}, CD8A={cd8a_exists}, CD8B={cd8b_exists}")
    
    cd4_mask = np.zeros(adata.shape[0], dtype=bool)
    cd8_mask = np.zeros(adata.shape[0], dtype=bool)
    
    if cd4_exists:
        cd4_expr = adata[:, 'CD4'].X.toarray().flatten()
        cd4_mask = cd4_expr > 1
        print(f"  CD4+细胞 (表达>1): {cd4_mask.sum()}个")
    
    if cd8a_exists or cd8b_exists:
        cd8_expr = np.zeros(adata.shape[0])
        if cd8a_exists:
            cd8_expr = np.maximum(cd8_expr, adata[:, 'CD8A'].X.toarray().flatten())
        if cd8b_exists:
            cd8_expr = np.maximum(cd8_expr, adata[:, 'CD8B'].X.toarray().flatten())
        cd8_mask = cd8_expr > 0.7
        print(f"  CD8+细胞 (表达>1): {cd8_mask.sum()}个")
    
    print(f"\n【保存barcodes文件】")
    
    all_barcodes = adata.obs_names.tolist()
    pd.Series(all_barcodes).to_csv(
        os.path.join(output_dir, f"{group_name}_cell_barcodes.txt"), 
        index=False, header=False
    )
    print(f"  ✓ 保存所有细胞barcode: {len(all_barcodes)}个")
    
    if cd4_mask.sum() > 0:
        cd4_barcodes = adata.obs_names[cd4_mask].tolist()
        
        cd4_df = pd.DataFrame({
            'index': range(len(cd4_barcodes)),
            'barcode': cd4_barcodes
        })
        cd4_df.to_csv(os.path.join(output_dir, f"{group_name}_CD4_barcodes.txt"), index=False)
        print(f" 保存CD4+细胞barcodes: {len(cd4_barcodes)}个")
        
    if cd8_mask.sum() > 0:
        
        cd8_barcodes = adata.obs_names[cd8_mask].tolist()
        
        cd8_df = pd.DataFrame({
            'index': range(len(cd8_barcodes)),
            'barcode': cd8_barcodes
        })
        
        cd8_df.to_csv(os.path.join(output_dir, f"{group_name}_CD8_barcodes.txt"), index=False)

for stage in stages:
    
    stage_cells = [cell for cell in hvg_expression.columns if cell_stages[cell] == stage]
    
    if stage_cells:
        # 筛选该阶段的数据
        stage_data = hvg_expression[stage_cells]
        
        # 保存为CSV
        output_file = path + f'{stage}_CD8_2000hvg.csv'
        stage_data.to_csv(output_file)
        
        print(f"  细胞数: {len(stage_cells)}个")
        print(f"  保存路径: {output_file}")
