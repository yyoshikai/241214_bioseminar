# 2 Materials and Methods
## 2.1 Dataset Curation
The datasets used to train RosettaFold All-Atom can be broadly grouped into three categories: protein-only, nucleic acid, and small molecule datasets. In this subsection, we describe the curation of each dataset. During training, protein sequence clusters are sampled within each dataset (further discussion of how often each dataset is sampled in subsection 2.6). Multiple Sequence Alignments (MSAs) and templates were generated as in RoseTTAFold2 (RF2) [11]. For MSAs, hhblits [75] was ran at successive e-value cutoffs until 10000 unique sequences with >50% coverage were found. For templates, HHsearch [76] was ran to find a maximum of 500 templates with probability >5%.
    RosettaFold All-Atomの学習に使用したデータセットは、タンパク質のみ、核酸、低分子の3つのカテゴリーに大別できる。このサブセクションでは、各データセットのキュレーションについて説明する。トレーニング中、タンパク質配列クラスターは各データセット内でサンプリングされる（各データセットのサンプリング頻度については2.6節で詳しく説明する）。MSAとテンプレートはRoseTTAFold2(RF2)[11]と同様に作成した。MSAについては、hhblits [75]を、50%以上のカバレッジを持つユニークな配列が10000個見つかるまで、連続したe-valueカットオフで実行した。テンプレートについては、HHsearch [76]を実行し、最大500個のテンプレートが確率5%以上で見つかるようにした。

### 2.1.1 Protein-Only Datasets
Datasets with protein-only examples come from both the Protein Data Bank (PDB) [77] and the AlphaFold2 (AF2) [1] distillation set.
    タンパク質のみの例を含むデータセットは、Protein Data Bank (PDB) [77]とAlphaFold2 (AF2) [1]の蒸留セットの両方から得られます。
#### PDB Monomers and Protein Complexes
Similar to RF2 [11], we train on protein monomer and protein complexes structures deposited into the PDB before April 30, 2020 with resolution below 4.5Å. For each chain, we find all contacting chains in all bioassemblies and featurize pairs of homomeric and heteromeric proteins by cropping around the interface. In cases where heteromeric complexes are from the same organism, we provide paired MSAs.
    RF2 [11]と同様に、2020年4月30日以前にPDBに寄託された4.5Å以下の解像度を持つタンパク質モノマーとタンパク質複合体の構造で学習します。 各鎖について、すべてのバイオアセンブリで接触しているすべての鎖を見つけ、界面周辺を切り取ることによって、ホモマーおよびヘテロマータンパク質のペアを特徴付けます。 ヘテロマー複合体が同じ生物由来の場合、ペアのMSAを提供します。

#### AlphaFold 2 Distillation Data 
We train RFAA on a set of UniRef50 structures predicted by AlphaFold2 in [78]. We follow RF2 and augment our training data with structures with AF2 mean pLDDT >70. We apply backbone losses for all structures but only apply sidechain losses for residues with AF2 pLDDT >90. The structures were predicted with only MSAs and no templates and we do not use templates for these examples in training.
    78]のAlphaFold2によって予測されたUniRef50の構造セットでRFAAのトレーニングを行う。RF2に従い、AF2平均pLDDT>70の構造でトレーニングデータを増強する。すべての構造に対してバックボーンロスを適用するが、AF2 pLDDT >90の残基に対してのみサイドチェーンロスを適用する。これらの構造はMSAのみで予測され、テンプレートは使用されていません。

### 2.1.2 Nucleic Acid Datasets
核酸データセット

We follow RoseTTAFold nucleic acid (RF-NA) [10] train on protein-nucleic acid complexes and RNA structures. The training setup is identical to the training set for RF-NA.
    RoseTTAFold核酸（RF-NA）[10]に従い、タンパク質と核酸の複合体やRNA構造を学習する。学習セットアップはRF-NAの学習セットと同じである。
### 2.1.3 Small Molecule Datasets
低分子データセット

This subsection deals with data that involves non-polymer small molecules (or non-linear polymers such as sugars), generally bound to a protein. Each item in the various small molecule datasets is represented by a selected subset of chains from a PDB entry. We build the dataset by iterating through each entry in the PDB with a bound, non-polymer chain. For each non-polymer chain, we firstdiscardthatentryif it is considered“non-biological”(seesubsection2.1.3).Otherwise, we treat that chain as a query ligand of interest and compute a list of every chain in that entry that contacts it (see subsection 2.1.3). We maintain a running list of transformations operate on each chain’s coordinates to place it into a global reference frame: this is necessary for symmetric assemblies in the PDB where a single chain may represent multiple copies of identical molecules in different,symmetricp o sitions.In allcases,hydrogensarenotm o deled.Finally,we compute MSAs and templates for each protein chain contacting each query ligand in the same manner as in subsection 2.1.1.
    このサブセクションは、一般的にタンパク質に結合した非高分子低分子（または糖のような非直線高分子）のデータを扱う。様々な低分子データセットの各項目は、PDBエントリから選択された鎖のサブセットで表される。データセットは、PDBの各エントリーから、結合している非ポリマー鎖を繰り返し検索することで構築される。各非高分子鎖について、まずそのエントリーが 「non-biological」(2.1.3節参照)であれば、そのエントリーを破棄します。そうでなければ、その鎖を目的のクエリーリガンドとして扱い、そのエントリーに含まれる、その鎖と接触するすべての鎖のリストを計算します(2.1.3節参照)。各鎖の座標をグローバルな参照フレームに配置するために、各鎖の座標に作用する変換の実行中のリストを保持します。これは、PDBの対称的な集合体では、1つの鎖が異なる対称的な位置にある同一分子の複数のコピーを表すことがあるため必要です。最後に、2.1.1節と同じ方法で、各クエリーリガンドに接触する各タンパク質鎖のMSAとテンプレートを計算します。

This definition of dataset item implies that a single PDB entry can comprise of multiple items in the dataset corresponding to different query l igands. For example, the PDB entry ‘5nag’ corresponds to two entries in our dataset: one for each bound ligand (FAD and 8R5). At training time, we crop each entry around the query ligand (further description in subsection 2.3.3), resulting in different contexts for the two query ligands.
    このデータセット項目の定義は、1つのPDB項目が、異なるクエリーリガンドに対応するデータセット内の複数の項目から構成されうることを意味している。例えば、PDBエントリ'5nag'はデータセット中の2つのエントリに対応する。トレーニング時に、クエリーリガンドを中心に各エントリーをクロップし（詳細は2.3.3節で説明）、その結果、2つのクエリーリガンドに対して異なるコンテキストが得られる。

We also note here that we only include in our small molecule datasets those PDB entries that are present in either the PDBBind [79] or BioLip [80] datasets.
    また、PDBBind[79]またはBioLip[80]データセットに存在するPDBエントリーのみを低分子データセットに含めることに注意する。
#### Ligand Filtering
リガンドのフィルタリング

We find that many of the non-polymer entities in PDB entries are non-biological: that is, they represent solvents or crystallization additives rather than binding partners found in a biological context. We curated a list of 3-letter component identifiers ( which have up to 3 letters) from the PDB corresponding to such non-biological molecules based on the BioLip database [80] and our own manual inspection as shown in Table 1. For example, GOL is the three letter code for glycerol, a common solvent.
    すなわち、それらは生物学的な文脈で見られる結合パートナーではなく、溶媒や結晶化添加物を表している。表1に示すように、BioLipデータベース[80]と我々の手作業による検査に基づいて、PDBからそのような非生物学的分子に対応する3文字の成分識別子（最大3文字）のリストを作成した。例えば、GOLは一般的な溶媒であるグリセロールの3文字コードです。
We remove all non-biological from our training set as they generally represent non-specific binding partners and/or molecules held in place by a crystal lattice rather than by protein interactions. We note that the list in Table 1 may not be exhaustive, but it filters out the most common examples we observed in the PDB.
    非生物学的分子は、一般的に非特異的結合パートナーや、タンパク質相互作用ではなく結晶格子によって保持されている分子を表すため、トレーニングセットからすべて削除した。表1のリストは完全ではないが、PDBで観察された最も一般的な例を除外している。



#### Contacting Chains
A protein chain is considered “in contact” with the query ligand if at least 5 Cα atoms of that protein chain are within 30Å of any atom in the query ligand or if any atom in the protein chain is within 5Å of any atom in the query ligand. This definition o f contacting protein chain allows for a broad biological context in which to predict the query ligand of interest.
    タンパク質鎖の少なくとも5個のCα原子がクエリーリガンドのどの原子からも30Å以内にある場合、あるいはタンパク質鎖のどの原子もクエリーリガンドのどの原子からも5Å以内にある場合、そのタンパク質鎖はクエリーリガンドと「接触」しているとみなされる。 このようなタンパク質鎖の定義により、クエリーリガ ンドを予測するための生物学的な背景が広くなる。
 
A non-polymer chain is considered “in contact” with the query ligand if at least one atom of that chain is within 5Å of the query ligand, or if all atoms are within 30Å of the query ligand. The model thus learns to predict not only the relative positions of a ligand and its protein binding partners, but also associated cofactors in the same or nearby binding pockets.
    非ポリマー鎖は、その鎖の少なくとも1つの原子がクエリ リガンドから5Å以内にある場合、あるいはすべての原子がクエリ リガンドから30Å以内にある場合、クエリリガンドと「接触」しているとみなされる。このようにしてモデルは、リガンドとそのタンパク質結合パートナーの相対的な位置だけでなく、同じ結合部位やその近傍の結合部位にある関連する補酵素も予測するように学習する。

#### Covalent Modification Dataset 
We separate query ligands that are covalently bonded to a protein into a distinct and separate dataset. For each such query ligand, we filtered out those that involved a covalent bond between an oxygen atom on the ligand and an oxygen atom on the protein, as well as protein-ligand fluorine-flourine bonds, covalent bonds to hydrogen atoms, and all cases where the length of the protein-ligand covalent bond was less than 1Å. We find that such covalent protein-ligand bonds are usually between a protein and a non-biological small molecule, as described in subsection 2.1.3.
    タンパク質と共有結合しているクエリーリガンドを別個のデータセットに分離した。各クエリーリガンドについて、リガンド上の酸素原子とタンパク質上の酸素原子の共有結合、タンパク質-リガンド間のフッ素-フッ素結合、水素原子への共有結合、タンパク質-リガンド間の共有結合の長さが1Å未満のものなどを除外した。 このようなタンパク質-リガンドの共有結合は、2.1.3節で述べたように、通常、タンパク質と非生物学的低分子との間の結合であることがわかった。

#### Metal Ion Dataset
We also separate query atoms that represent metal ions into a distinct dataset. In order to avoid expanding the vocabulary set of our model too greatly and potentially introducing non-specific metal binding sites into the dataset, we only train on the metals whose 3-letter PDB codes are listed in Table 2.
    また、金属イオンを表すクエリー原子を別個のデータセットに分離する。モデルの語彙を増やしすぎたり、非特異的な金属結合部位をデータセットに取り込む可能性を避けるため、PDBコードの3文字が表2に記載されている金属のみを対象とした。

#### Small Molecule Datasets
Every remaining query ligand can be grouped into one of the fol- lowing groups: multi-residue, multi-protein assembly, and single-protein assembly. We describe the three sub-datasets here:
    残りの全てのクエリーリガンドは、以下のグループの1つにグループ化できる：多残基、多タンパク質アセンブリー、単一タンパク質アセンブリー。ここでは3つのサブデータセットについて説明する：
1. Protein/Small Molecule Complex: Dataset containing all single residue small molecules that only interact with a single other protein chain. However, the query ligand in each item of this dataset lies on a single residue in the PDB entry.
    1. タンパク質/低分子複合体： 他の1つのタンパク質鎖としか相互作用しない1残基の低分子をすべて含むデータセット。ただし、このデータセットの各アイテムのクエリーリガンドは、PDBエントリーの1残基に存在する。
2. Protein/Multi-Residue Ligand Complex: Dataset containing all “multi-residue” lig- ands. Such ligands exist as multiple residues (or a single residue with multiple bioassembly transforms) in their respective cif file, and usually represent sugar chains or small peptides.
    2. タンパク質/多残基リガンド複合体： すべての 「multi-residue 」リガンドを含むデータセット。このようなリガンドは、それぞれのcifファイル中で複数の残基（または複数の生体組織化変換を伴う単一残基）として存在し、通常は糖鎖や小さなペプチドを表す。
3. Protein/Small Molecule Assembly: Dataset containing non-protein biomolecular con- text (small molecules, covalent modifications, multi-residue ligands, metal ions) and > 1 protein chain. This dataset is distinct from 1. because we generate paired MSAs when there are two or more distinct (heteromeric) protein chains in a complex.
    3. タンパク質/低分子アセンブリ： タンパク質以外の生体分子（低分子、共有結合修飾、多残基リガンド、金属イオン）と1本以上のタンパク質鎖を含むデータセット。このデータセットは1.とは異なり、複合体中に2つ以上の異なる（ヘテロ）タンパク質鎖が存在する場合に、ペアMSAを生成する。

#### Dataset Clustering
We cluster each entry in each dataset by their primary protein partner. For protein-only datasets, this is either the protein monomer, or in the case of hetero-oligomers, we arbitrarily select the firstchainthatappearsinthePDBentryastheprimaryc hain.For nucleic acid and small molecule datasets, we designate as the primary protein partner the protein chain with the greatest number of atoms within 5Å of the query nucleic acid/ligand. All items in each database are then clustered using MMSeqs2 [81] using the default hyper-parameters and subsequently used for dataset sampling at training time.
    各データセットの各エントリーを主要なタンパク質パートナーによってクラスタリングした。核酸と低分子のデータセットでは、クエリー核酸/リガンドから5Å以内にある原子の数が最も多いタンパク質鎖を一次タンパク質パートナーとする。 各データベースの全アイテムは、MMSeqs2 [81]を用いてデフォルトのハイパーパラメータでクラスタリングされ、その後、トレーニング時のデータセットサンプリングに使用される。

#### Final Statistics
After all the filters w ere a pplied, w e k ept c hains t hat h ave r esolution < 4.5Å and that were deposited in the PDB before April, 30th, 2020 for the purposes of held-out evaluation on PDB entries from 2021 onward. 10% of protein sequence clusters were held out of training for validation during training. The final number of items and clusters in each dataset is shown in Table 3. The number of validation items and clusters in each dataset is shown in Table 4.
    すべてのフィルターを適用した後、解像度が4.5Å未満で、2020年4月30日以前にPDBに登録された配列を抽出し、2021年以降のPDBエントリーの評価用に保留した。 10%のタンパク質配列クラスターは、トレーニング中に検証のためにトレーニングから除外された。 各データセットの最終的なアイテム数とクラスター数を表3に示す。 各データセットの検証項目数とクラスタ数を表4に示す。

### 2.1.4 CSD Dataset
2.1.4 CSDデータセット

In addition to molecules in the PDB, we augment our training dataset with small molecule crystal structures from the Cambridge Structure Database (CSD v5.43; November 2021) [12]. We filter structures based on the following metrics: 1) resolution < 5Å, 2) not polymeric, 3) greater than 5 atoms resolved, 4) less than 100 atoms resolved, and 5) ability to be parsed by OpenBabel [82]. We sample molecules with equal probability and separate a validation set that does not have Tanimoto score >0.75 to any molecule in training.
    PDBの分子に加えて、Cambridge Structure Database (CSD v5.43; November 2021) [12]の低分子結晶構造をトレーニングデータセットに追加します。以下のメトリクスに基づいて構造をフィルタリングする： 1) 分解能が5Å未満、2) 重合していない、3) 分解能が5原子以上、4) 分解能が100原子未満、5) OpenBabel [82]で解析可能。 等確率で分子をサンプリングし、Tanimoto スコアが0.75 を超えない検証セットをトレーニングのどの分子とも分離する。

### 2.1.5 Negative Datasets
負のデータセット

Following RF2 and RF2NA, we use a set of “negative” interactions to help the network focus on relevant features that constitute binding. Briefly, for proteins this means showing examples where we randomly pair chains and only assess a loss on each individual chain and for nucleic acids it involves mutating bases that make essential contacts for the formation of the complex. Negative examples were not shown for any small molecule dataset.
    RF2とRF2NAに続いて、結合を構成する関連する特徴にネットワークが集中できるように、「負の 」相互作用のセットを使う。簡単に言うと、タンパク質の場合は、鎖をランダムにペアにして、個々の鎖の損失だけを評価する例を示すことを意味し、核酸の場合は、複合体の形成に不可欠な接点を作る塩基を変異させることを含む。ネガティブな例は、どの低分子データセットでも示されなかった。

## 2.2 Modeling Arbitrary Biological Inputs
Architectures for modeling both protein structures and nucleic acid structures have been previously described. In this subsection, we describe the necessary changes to such architectures to model small molecules, covalent modifications to protein structures and arbitrary non-canonical amino acids.
    タンパク質構造と核酸構造の両方をモデリングするためのアーキテクチャーは、以前に説明されている。この小節では、低分子、タンパク質構造への共有結合修飾、および任意の非正規アミノ酸をモデル化するために、そのようなアーキテクチャに必要な変更について説明する。

### 2.2.1 Expanded Input Sets

The most significant architectural change from existing protein structural networks is the expanded input features that RFAA takes in, which we describe in the following subsections.
    既存のタンパク質構造ネットワークとの最も大きな違いは、RFAAの入力機能の拡張である。

#### New Tokens
The original RoseTTAFold architecture had 22 tokens: 20 amino acids, 1 unknown and 1 gap token. The RF-Nucleic Acid expanded this token set by 10, adding 8 distinct tokens for the 4 DNA and 4 RNA bases, and 2 tokens for unknown DNA base and unknown RNA base, respectively. The RFAA architecture includes 46 additional tokens representing individual atoms with the element types shown in Table 5, a token for deprotonated histidine (unused in practice, left in for legacy reasons) and an unknown atom token for a total token count of 80.
    オリジナルのRoseTTAFoldアーキテクチャには22のトークンがあった：20個のアミノ酸、1個の未知、1個のギャップトークン。RF-Nucleic Acidはこのトークンセットを10個拡張し、4個のDNA塩基と4個のRNA塩基を表す8個のトークンと、それぞれ未知のDNA塩基と未知のRNA塩基を表す2個のトークンを追加した。RFAAのアーキテクチャには、表5に示す要素タイプを持つ個々の原子を表す46のトークンが追加され、脱プロトン化ヒスチジンのトークン（実際には未使用、レガシー上の理由で残っている）、未知の原子のトークンが加わり、合計トークン数は80となった。

#### Bond Connectivity
When predicting structures of arbitrary molecules, it is important for the network to know the bond connectivity of those molecules. To provide this information to the network, we pass in the bond connectivity of input molecules as a 2D bond adjacency matrix as an input. We designate 7 bond types representing single bond, double bond, triple bond, aromatic bond, residue-residue (or base-base), residue-ligand atom bond and other bond type. In practice, the “other” bond type is not used but exists for historical reasons. The residue-residue “bond” type exists to be able to provide bond features for protein and nucleic acid inputs, so that the input bond matrix is always of the same dimension as the input. Residue-atom bonds exist in order to do a process we call residue atomization, which is used to model ligands that are covalently bonded to a residue and arbitrary non-canonical amino acids. This process is described further in subsection 2.3.4.
    任意の分子の構造を予測する場合、ネットワークはそれらの分子の結合結合性を知ることが重要である。この情報をネットワークに提供するために、入力分子の結合結合性を2次元結合隣接行列として入力として渡す。単結合、二重結合、三重結合、芳香族結合、残基-残基（または塩基-塩基）結合、残基-リガンド原子結合、その他の結合の7種類を指定する。実際には「その他」の結合型は使用されないが、歴史的な理由から存在する。残基-残基結合タイプは、入力結合行列が常に入力と同じ次元になるように、タンパク質や核酸の入力に対して結合の特徴を提供できるようにするために存在する。残基-原子結合は、残基と任意の非共有結合アミノ酸に共有結合しているリガンドをモデル化するために使用される残基原子化と呼ぶプロセスを実行するために存在する。この処理については、2.3.4節で詳しく説明する。

#### Chiral Features
Another key bias in more generalized biomolecular modeling is chirality. Afore- mentioned features like atom types and bond connectivities are not sufficient to specify the chi- rality of input molecules, so we provide chiral features explicitly to the network at each chiral center. In this work, we only deal with tetrahedral chirality and leave more complicated forms of stereochemistry to future work.
    より一般化された生体分子モデリングにおけるもう一つの重要なバイアスはキラリティーである。前述した原子の種類や結合の連結性のような特徴は、入力分子のキラリティーを指定するには不十分であるため、我々は各キラル中心でネットワークに明示的にキラリティーの特徴を与える。本研究では、四面体のキラリティのみを扱い、より複雑な立体化学は今後の研究に委ねる。

For each tetrahedral chiral center, we enumerate all sets of three heavy atom neighbors in all orders. For each ordering, we compute the pseudo-dihedral angle between those four points (center and 3 heavy atoms) and note whether that angle is positive or negative, which determines the chirality of the center uniquely. For ideal tetrahedral geometry, the magnitude of the dihedral angle is √arcsin(1/ 3). See subsection 2.4.5 for further details.
    各四面体キラル中心について、3つの重原子の隣接集合をすべての次数で列挙する。各順序について、これらの4点（中心と3つの重原子）間の擬二面角度を計算し、その角度が正か負かを記録する。理想的な四面体形状の場合、二面角の大きさは√arcsin(1/3)である。詳細は2.4.5節を参照。

We compute the difference b etween t he p redicted d ihedral a ngle a nd t he i deal d ihedral a ngle of the chiral center, and pass the coordinate-wise gradients of the difference a s a n i nput t o t he 3D track of the network. This representation of chirality provides direct signal to the network on how to update atomic coordinate positions in order to obey the ideal tretrahedral geometry. This process is described further in subsection 2.4.5.
    我々は、キラル中心の再描画された二面体アングルと取引された二面体アングルとの差を計算し、その差の座標ごとの勾配をネットワークの3Dトラックに入力として渡す。このキラリティの表現は、理想的な三面体幾何学に従うために、原子座標位置をどのように更新するかという直接的なシグナルをネットワークに与える。このプロセスについては2.4.5節で詳しく説明する。

#### Atom Frames
Key to the success of AF2 and RF2 is the frame aligned point error loss [1]. This involves aligning N-Cα-C backbone frames of predicted structures to true structures and then measuring the error of all the other predicted atoms with respect to the true structure in that alignment. This loss has attractive properties for biomolecules such as not being invariant to reflections w hich a llows t he n etwork t o p redict c orrect c hirality. W e c onstruct c anonical frames for each atom in small molecules comprising of atoms and their bonded neighbors. We achieve this by iterating through all bonded triplets of atoms and assigning each triplet a priority based on the bonded atoms, depicted in Table 7. The process for constructing canonical frames from a ligand is outlined below:
    AF2とRF2の成功の鍵は、フレーム整列点誤差損失である[1]。これは、予測された構造のN-Cα-Cバックボーン・フレームを真の構造に整列させ、その整列における他のすべての予測原子の真の構造に対する誤差を測定するものである。この損失は生体分子にとって魅力的な特性を持っており、例えば反射に対して不変でないため、ネットワークが正しいヒラリティを再認識することができる。我々は、原子とその結合した近傍原子からなる低分子の各原子について、非対称フレームを構築する。これは、表7に示すように、原子の結合した三重項をすべて繰り返し、各三重項に結合原子に基づく優先順位を割り当てることで実現する。リガンドからカノニカルフレームを構築するプロセスを以下に概説する：


1. Construct a graph where each node is an atom and each edge is a bond.
    1. 各ノードが原子で、各辺が結合であるグラフを作る。

2. For each atom in the graph:
    2. グラフの各原子について
    (a) Enumerate through all paths of length three containing that atom. If there exist paths such that the given atom is in the center, exclude all other paths.
        (a) その原子を含む長さ3のすべてのパスを列挙する。与えられた原子が中心にあるようなパスが存在する場合、他のパスはすべて除外する。
    (b) Compute frame priorities for each atom in each path and make a list of frame priorities in increasing order.
        (b) 各パスの各原子についてフレームの優先度を計算し、フレームの優先度を昇順に並べたリストを作成する。
    (c) For each such path of length three, sort them by lexicographic atom frame priority, e.g. for two paths A and B, path A will appear before path B if and only if either the lowest frame priority in path A is less than the lowest frame priority in path B, or they are equal and the second lowest frame priority in path A is less than the second lowest in path B, and so on.
        (c)このような長さ3の各パスについて、アトムのフレーム優先度を辞書順に並べ替える。例えば、2つのパスAとBについて、パスAのフレーム優先度がパスBのフレーム優先度の最低値より小さいか、あるいはそれらが等しく、かつパスAのフレーム優先度の2番目に低いものがパスBのフレーム優先度の2番目に低いものより小さい場合に限り、パスAはパスBより先に現れる。
    (d) The first path in the lexicographic order is chosen as the frame for this atom, and the order of the frame is determined by frame priority in increasing order.
        (d)辞書順の最初のパスがこの原子のフレームとして選択され、フレームの順序はフレームの優先順位が高い順に決定される。

This process deterministically computes a local coordinate frame for each atom in arbitrary molecules (with at least 3 atoms). If a frame has an unresolved atom, it still is assigned as the canonical frame but is not used in loss calculation. The usage of the atom frames is further discussed in subsections 2.4.5, 2.5.5 and 2.5.4. Importantly, these features are constructed from the bond graph so they maintain the permutation invariance of the inputs to the network.
    このプロセスは、任意の分子（少なくとも3つの原子を持つ）の各原子のローカル座標フレームを決定論的に計算する。フレームに未解決の原子がある場合、そのフレームは正準フレームとして割り当てられますが、損失計算には使用されません。原子フレームの使用法については2.4.5、2.5.5、2.5.4節で説明する。重要なことは、これらの特徴は結合グラフから構築されるため、ネットワークへの入力の順列不変性を維持することである。

#### Positional Encodings
To break permutation symmetry for sequences, we use a signed relative positional encoding for protein sequences and nucleic acid bases [1]. Atomic graphs require permutation symmetry so we do not provide a relative positional encoding. For atomic inputs, we provide a separate embedding that measures the shortest distance between any pair of atoms in the bond graph described in subsection 2.2.1. We develop a generalization of these two embeddings for cases where atom nodes are bonded to residues to encode the distance between an atom and its closest bonded residue. Further details are provided in subsection 2.4.1.
    配列の順列対称性を破るために、タンパク質配列と核酸塩基の符号付き相対位置エンコーディングを用いる[1]。原子グラフは順列対称性が必要なので、相対位置エンコーディングは提供しない。原子入力に対しては、2.2.1節で説明した結合グラフの任意の原子のペア間の最短距離を測定する埋め込みを別に提供する。原子ノードが残基と結合している場合のために、これら2つの埋め込みを一般化し、原子と最も近い結合残基間の距離をエンコードする。詳細は2.4.1節で述べる。

## 2.3 Data Pipeline
Our data pipeline involves taking raw data from cif files and formatting the data into input tensors for the network. We first will describe the inputs to the network and then go through details of how each dataset is preprocessed. We use OpenBabel to parse the ideal sdf files for each PDB but do not use any chemical quantities computed by it (just element types, bond types and whether an atom is chiral center). We find this to be preferable because OpenBabel can process all ideal sdf files provided by the PDB so we do not exclude examples because of parsing e rrors. We compute the direction of the chiral center as discussed in subsection 2.4.5.
    我々のデータパイプラインは、cifファイルから生データを取り出し、データをネットワークの入力テンソルにフォーマットする。最初にネットワークへの入力を説明し、次に各データセットがどのように前処理されるかの詳細を説明する。各PDBの理想的なsdfファイルを解析するためにOpenBabelを使用するが、それによって計算された化学量は使用しない（元素タイプ、結合タイプ、原子がキラル中心であるかどうかだけ）。OpenBabelはPDBから提供されたすべての理想的なsdfファイルを処理できるため、パースエラーによって例を除外することがないためです。2.4.5節で説明したように、キラル中心の方向を計算する。

### 2.3.1 Inputs for RFAA
Remaining features such as MSAs and templates are handled identically for proteins to RF2. The coordinate dimension, 36, reflects t he m aximum a mount o f h eavy a toms a nd h ydrogens possible in a residue or base. The small molecule tokens are appended to the first s equence i n a ll the MSA features and the remaining MSA sequences are initialized with gap tokens. Small molecules receive empty template features which are concatenated block diagonally to the protein features. A detailed description of the inputs are shown in Table 8.
    MSAやテンプレートなどの残りのフィーチャーは、RF2に対するタンパク質と同じように扱われます。座標寸法36は、残基や塩基に含まれる重金属や水素の最大数を反映する。低分子のトークンはすべてのMSAフィーチャーの最初の配列に付加され、残りのMSA配列はギャップトークンで初期化されます。低分子は空のテンプレート特徴量を受け取り、タンパク質特徴量に斜めに連結される。入力の詳細を表8に示す。

### 2.3.2 Featurization of Symmetric Permutations
When featurizing multimers (empirically in the homomer and small molecule assembly datasets), there are often identical chains that, if swapped, result in the identical complex. We want the gradient of the loss to push the network towards the relabeled complex that is closest to the pre- diction so during preprocessing we track which chains can be swapped so that we can deconvolute which permutation to apply the loss on during training. We use the same scheme to account for permutation swaps of atoms in small molecule structures.
    ホモマーデータセットや低分子アセンブリーデータセットで経験的に）多量体をフィーチャー化する場合、しばしば同じ鎖が存在し、それを入れ替えると同じ複合体になる。そのため、前処理において、どの鎖を入れ替えることができるかを追跡し、学習時にどの順列に損失を適用するかを決定できるようにする。我々は、低分子構造における原子の順列入れ替えを考慮するために、同じスキームを使用する。

### 2.3.3 Featurization of Protein Small Molecule Complexes
A similar preprocessing procedure was followed for the small molecule protein complex dataset, the multichain residue ligand dataset, the covalent modification dataset and the protein-small molecule assembly dataset (multiple contacting protein chains). Each training example centers around a single query molecule in a specific bioassembly. Based on the details of the bioassembly features for the nearest chains (both protein and other small molecules) are constructed. There are two types of biomolecular contexts that are sampled stochastically. First, metal ions in the presence of other small molecules are sampled stochastically because often is is not a priori known when solving a structure whether there will be a metal crystallized in the pocket. Second, if there modified residues present in the cif file, they are featurized as an atomized residues rather than their canonicalized version.
    低分子タンパク質複合体データセット、多鎖残基リガンドデータセット、共有結合修飾データセット、タンパク質-低分子アセンブリーデータセット（複数の接触タンパク質鎖）についても同様の前処理を行った。各トレーニング例は、特定のバイオアセンブリーにおける単一のクエリー分子を中心としている。バイオアセンブリーの詳細に基づいて、最も近い鎖（タンパク質と他の低分子の両方）の特徴が構築される。確率的にサンプリングされる生体分子コンテキストは2種類ある。第一に、他の低分子の存在下にある金属イオンは確率的にサンプリングされる。なぜなら、構造を解くときに、ポケットに結晶化した金属が存在するかどうかが事前にわからないことが多いからである。第二に、cifファイルに修飾残基がある場合、それらは正準化されたバージョンではなく、原子化された残基としてフィーチャーされます。

Due to memory restrictions, we then perform a cropping procedure to select a subset of nodes to represent a training example. The cropping procedure samples a random atom on the query molecule and computes the distance to all other atoms or Cα atoms in proteins. It then selects the top n_crop nodes to include in the crop. In our early experiments we found that sometimes, there were protein chains that were either too short or far away from the ligand in euclidean space that were included in the crop with insufficient context to be predicted accurately. In these cases, the gradient was dominated by the incorrect prediction of those protein chains and not on the correct docking of the small molecule. To remedy this, after finding the top n_crop nodes, we iterate through all the protein chains that were in the crop and remove any chains with <10 contacts to other nodes in the crop or with less than 8 residues. After these chains are removed, we noticed that certain molecules in the crop also did not have sufficient contacts to be docked so we iterate through all the molecules in the crop and remove any molecules that have <4 contacts to a protein chain. The exact logic is shown in 1.
    メモリの制約から、次に、学習例を表すノードのサブセットを選択するために、切り出し手順を実行する。クロッピング手順は、クエリ分子上のランダムな原子をサンプリングし、タンパク質中の他のすべての原子またはCα原子との距離を計算する。そして、上位n_cropノードを選択してクロッピングに含める。初期の実験では、ユークリッド空間においてリガンドから短すぎたり遠すぎたりするタンパク質鎖が、正確な予測には不十分な文脈でクロップに含まれることがあることがわかった。このような場合、勾配は、低分子の正しいドッキングではなく、それらのタンパク質鎖の不正確な予測に支配されていた。これを改善するために、上位n_cropノードを見つけ た後、クロップ内にあったすべてのタンパク質鎖を繰り返し検 索し、クロップ内の他のノードとの接触が10未満、または残基が8 未満の鎖を削除した。これらの鎖が除去された後、我々は作物中の特定の分子もドッキングするのに十分な接触を持っていないことに気づいたので、作物中のすべての分子を繰り返し処理し、タンパク質鎖との接触が4未満の分子をすべて除去する。正確なロジックを1に示す。

It is evident from this cropping process that full subunits in large symmetric assemblies could be cropped out. Since we compute potential symmetric relabeling of chains before cropping (see sub- section 2.3.2), certain chains that were computed as potential symmetric relabelings are no longer valid (specifically because the small molecule context should drive the network to predict a spe- cific interface when a symmetric oligomer could have multiple distinct protein-protein interfaces). After cropping, we reiterate through the precomputed symmetric permutations and remove those that are no longer possible given the chains that were removed during cropping.
    この切り出しプロセスから、大きな対称的集合体の完全なサブユニットが切り出される可能性があることは明らかである。我々はクロッピングの前に鎖の潜在的な対称的な再配置を計算するので（2.3.2節参照）、潜在的な対称的な再配置として計算されたある鎖は、もはや無効である（特に、対称的なオリゴマーが複数の異なるタンパク質-タンパク質界面を持つ可能性がある場合、低分子のコンテキストは、特定の界面を予測するようにネットワークを駆動する必要があるため）。クロッピングの後、事前に計算された対称順列を再計算し、クロッピング中に削除された鎖を考慮して、もはや不可能なものを削除する。

### 2.3.4 Featurization of Atomized Protein Examples
Training examples for atomized proteins first are featurized identically to protein monomer exam- ples (except the stochastic homomer featurization is turned off, see subsection 2.6). After cropping, a number of residues is sampled from Uniform(3,5) and that number of (fully resolved) contiguous residues is chosen for atomization. If there are not enough valid residues in the crop to atomize, we treat the example as a monomer example. We then take that selection of residues, featurize them using all the small molecule features (atom tokens, bond features, chirality inputs). We also provide bond tokens to indicate bonds between the first N token in the atomized region to the previous residue and the last C token to the following residue. Finally, the MSA and template information for these residues is removed from the input features so the network must learn how to generate their structures and poses from the atomic information. Practically, we precompute the atoms, bonds and chiralities of each atom in each residue and convert the features as shown in Algorithm 2. Symmetric swaps of sidechain atoms are accounted for in the same manner as subsection 2.3.2.
    アトマイズタンパク質のトレーニング例は、まずタンパク質モノマーと同じようにフィーチャライズされます（ただし、確率的ホモマーフィーチャライズはオフにされます、サブセクション2.6を参照）。切り出し後、残基の数がUniform(3,5)からサンプリングされ、その数の（完全に分解された）連続残基が原子化に選ばれます。アトマイズするのに十分な残基がクロップ内にない場合、その例をモノマーの例として扱います。次に、選択した残基を、すべての低分子特徴（原子トークン、結合特徴、キラリティ入力）を使用してフィーチャライズします。また、原子化領域の最初のNトークンから前の残基、最後のCトークンから次の残基までの結合を示す結合トークンを提供します。最後に、これらの残基のMSAとテンプレート情報は入力特徴から削除されるため、ネットワークは原子情報から残基の構造とポーズを生成する方法を学習する必要があります。実際には、各残基の原子、結合、キラリティーを事前に計算し、アルゴリズム2に示すように特徴量を変換する。側鎖原子の対称的な入れ替わりは2.3.2節と同じ方法で説明する。

### 2.3.5 Featurization of Covalently Bound Ligands
Covalently bound ligands are preprocessed very similarly to small molecule complexes. The cova- lent modifications are modelled just as other ligands would be. The residue that has the covalent bond to the ligand is atomized using the same method as subsection 2.3.4. The atom in the ligand and the atom in the atomized residue is provided in the bond features (eg. single bond between atom i from modification and atom j in atomized residue). All other featurization (protein MSA, templates etc) remains the same as other datasets.
    共有結合したリガンドは、低分子複合体と非常に同様に前処理される。共有結合の修飾は、他のリガンドと同様にモデリングされる。リガンドと共有結合している残基は、2.3.4節と同じ方法で原子化される。リガンドの原子と原子化残基の原子は、結合フィーチャーで提供されます（例：修飾の原子iと原子化残基の原子jの間の単結合）。他のすべての特徴化（タンパク質MSA、テンプレートなど）は、他のデータセットと同じです。

### 2.3.6 Featurization of Metal Ions
Metal Ions are provided to the network as a single atom ligand. The only difference is that since metal ions only have a single atom, they do not have their own canonical frame. In these cases, the network does not receive a frame input and there is no loss calculated with respect to the frame of the ion (there are still gradients from the error of the placement of the ion with respect to the other frames in the structure).
    金属イオンは単一原子のリガンドとしてネットワークに提供される。唯一の違いは、金属イオンは原子を1つしか持たないため、独自のカノニカルフレームを持たないことである。このような場合、ネットワークはフレーム入力を受け取らず、イオンのフレームに関して計算される損失はありません（構造内の他のフレームに対するイオンの配置の誤差による勾配は残っています）。

### 2.3.7 Featurization of CSD Small Molecule Crystals
Asymmetric units of crystal lattices from the CSD are featurized identically to small molecules that are bound to proteins. The network is then tasked with predicting the atomic coordinates of the molecule. Molecules with less than 5 atoms, greater than 100 atoms, polymers or resolution >5Å are discarded from the dataset. Remaining molecules that could be parsed by OpenBabel were used to train the network.
    CSDの結晶格子の非対称単位は、タンパク質に結合している低分子と同じように特徴づけられる。その後、ネットワークは分子の原子座標を予測する。原子数が5未満の分子、原子数が100を超える分子、ポリマー、または分解能が5Åを超える分子はデータセットから除外される。 OpenBabelで解析できた残りの分子をネットワークの学習に使用した。


































































