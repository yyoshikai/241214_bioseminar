要約
    AlphaFold2を開発した。
    類似した配列がないタンパク質についても, 原子スケールの精度での構造予測を達成した。
Introduction
    タンパク質の構造予測には, 物理化学的な相互作用の計算ベースの手法と, 進化の情報を使う手法がある。
    以前は前者が主流だったが, 計算量が多すぎてうまくいかなかった。
    近年は後者が進展している。
    しかしいずれの手法も実用的な精度にはなっていない。

    本研究では, 原子サイズレベルの精度で予測可能なモデルを初めて開発した。
    開発したモデルは, CASP14で予測精度に関する多くの指標で2位のモデルを大きく上回った。強力なテンプレートが利用可能な場合でも, テンプレートベースの手法を上回った。domain-packingを持つ長いポリペプチドにも適用できた。
    AlphaFold2の学習データ収集後にPDBBindに登録されたタンパク質についても正確に構造を予測しており, CASP14の結果がより大きく一般的なデータセットでも成り立つことを示した。

The AlphaFold Network
    AlphaFoldは, 進化的, 物理的, 幾何学的な制約に基づいた新しいニューラルネットワーク構造と学習手法により, 高いタンパク質構造予測の精度を達成した。
    特に,
        MSAとpairwise featureをまとめて埋め込む構造
        正確なend-to-endの構造予測を可能にする, 新しい出力の表現形式とそのロス関数
        同変性のある新しいattention機構
        予測を繰り返し精緻化するためのintermediate lossの使用
        構造とともに学習させるmasked MSAA loss
        モデルが自分て自分を蒸留し, 自分の精度を予測することによる正解構造のないアミノ酸配列からの学習
    を示す。

    AlphaFoldは, タンパク質のアミノ酸配列とそのホモログの配列から, タンパク質中の全ての重原子の座標を直接予測する。
        データベースや, (アミノ酸配列が与えられたときの?)MSAの構築方法についてはMethodを参照
        *回転はどうしているのか？
        以下ではモデルの最も重要なアイデアと構成要素を説明する。
        全てのネットワーク構造と学習手法についてはMethodを参照

    まず, Evoformerというブロック数層を通じ, 入力から$N_{seq}\times N_{res}$の表現($N_{seq}$は多分MSAの数)と, $N_{res}\times N_{res}$の残基間の相互作用の表現を生成する。
    `In`

# 原論文と解説

AlphaFold2は, backboneとhead的なものからなる。
backboneは, MSAのembeedding((B*?)N_msa*N_res*channel)と, 残基のペアの表現(N_res*N_res, 以下ペア表現)を入力とし, 同じものを出力とする。
backboneは, 48層のEvoformerからなり, それぞれが上記の入出力を受け取り出力する。
1層のEvoformerは, 
    (1)  MSAの表現をN_res方向に(つまり, 各MSAの表現について)self-attention。このときペア表現(を適当に変換したもの？)をバイアスとして加える。 この操作の前後はresidual connection付き
    (2) MSAの表現をN_msa方向にself-attention。residual connection付き
    (3) Transition...feedforwardのこと？
    (4) MSA表現の外積を取り(...どんな操作?), ペア表現に加える
    (5) ペア表現についてそれぞれ更新(図を参照...)
で構成される。
headは, MSA表現のうち大元のアミノ酸配列に対応するものと, pair表現を入力とする。共通するある一点からのそれぞれの残基の平行移動・回転がどの程度あるかでタンパク質構造を表現し, それを生成する。なので, 生成した構造がちゃんと鎖にならないこともあるが, それはfine tuningの段階で学習する。最終的にはAmber field?のエネルギー計算でつなげる(?)
headはよくわからなかった。同変性を持つように色々と工夫しているらしい。


## Evoformer
    本ネットワークの構成要素であるEvoformerの基本方針は, タンパク質構造予測を3次元空間上のグラフを推論する問題と捉えることである。
    ペア表現の要素は残基間の関係についての情報を表し, MSA表現の列(1残基の全MSAについての情報)は各残基の情報を, 行(1つのMSAの全残基の情報)はそれらが出現する配列の情報を表す。この枠組みの中で, いくつかの異なるアップデートを定義した。
    ペア表現は, MSA表現の外積(列方向に和をとる)によってアップデートされる。この処理は, これまでの研究では一度だけでなくブロックごとに毎度行われ, 更新されていくMSA表現からペア表現への継続的な情報伝達を可能にする。
The key principle of the building block of the network — named Evo-former (Figs. 1e, 3a) - is to view the prediction of protein structures as a graph inference problem in 3D space in which the edges of the graph are defined by residues in proximity. The elements of the pair representation encode information about the relation between the residues (Fig. 3b). The columns of the MSA representation encode the individual residues of the input sequence while the rows represent the sequences in which those residues appear. Within this framework, we define a number of update operations that are applied in each block in which the different update operations are applied in series.
The MSA representation updates the pair representation through an element-wise outer product that is summed over the MSA sequence dimension. In contrast to previous work , this operation is applied within every block rather than once in the network, which enables the continuous communication from the evolving MSA representation to the pair representation.
    Evo-former（図1e、図3a）と名付けられたこのネットワークの構成要素の主要な原理は、タンパク質構造の予測を3次元空間におけるグラフ推論問題として捉えることである。ペア表現の要素は、残基間の関係に関する情報をコード化する（図3b）。MSA表現の列は入力配列の個々の残基を表し、行はそれらの残基が現れる配列を表す。このフレームワークの中で、異なる更新操作が直列に適用される各ブロックで適用されるいくつかの更新操作を定義する。

    ペア表現内部では2種類のアップデートがあり, どちらも三角不等式などペア情報を単一の3次元構造として表すのに必要な制約から発想したものである。具体的には, グラフ上の全ての3頂点の組について, それがなす三角形の一辺の情報を残り2辺の情報でアップデートする。これは元々attention(図の後半2つ?)の軽量な代替物として考案されたが, 両方使った方が精度が上がったのでそうした。attentionはつまりペア表現についてのaxial attentionになると思う。
Within the pair representation, there are two different update patterns. Both are inspired by the necessity of consistency of the pair representation—for a pairwise description of amino acids to be representable as a single 3D structure, many constraints must be satisfied including the triangle inequality on distances. On the basis of this intuition, we arrange the update operations on the pair representation in terms of triangles of edges involving three different nodes (Fig. 3c). In particular, we add an extra logit bias to axial attention31 to include the ‘missing edge’ of the triangle and we define a non-attention update operation ‘triangle multiplicative update’ that uses two edges to update the missing third edge (see Supplementary Methods 1.6.5 for details). The triangle multipli- cative update was developed originally as a more symmetric and cheaper replacement for the attention, and networks that use only the attention or multiplicative update are both able to produce high-accuracy structures. However, the combination of the two updates is more accurate.
    ペア表現では、2つの異なる更新パターンがある。どちらもペア表現の一貫性の必要性から着想を得たもので、アミノ酸のペア表現を1つの3次元構造として表現するためには、距離に関する三角形の不等式を含む多くの制約を満たす必要がある。この直感に基づき、ペア表現の更新操作を3つの異なるノードを含むエッジの三角形で整理した（図3c）。特に、三角形の「欠けている辺」を含めるために、軸方向の注意31に余分なロジットバイアスを追加し、欠けている3番目の辺を更新するために2つの辺を使用する非注意の更新操作「三角形乗法更新」を定義する（詳細は補足方法1.6.5を参照）。三角形乗法更新は、もともと注意に代わる、より対称的で安価な更新として開発されたものであり、注意または乗法更新のみを使用するネットワークは、どちらも高精度な構造を生成することができる。しかし、2つの更新を組み合わせた方がより正確である。

    また, MSA表現内ではaxial attentionの変種として, ペア表現からのlogit(...バイアス?)をバイアスとして追加する。これにより, ペア表現とMSA表現を混合して処理することができる。
We also use a variant of axial attention within the MSA representation. During the per-sequence attention in the MSA, we project additional logits from the pair stack to bias the MSA attention. This closes the loop by providing information flow from the pair representation back into the MSA representation, ensuring that the overall Evoformer block is able to fully mix information between the pair and MSA representations and prepare for structure generation within the structure module.
    また、MSA 表現の中で軸方向注意の変種を使用する。MSAのシーケンスごとの注意の間に、ペアスタックから追加のロジットを投影し、MSAの注意を偏らせる。これにより、ペア表現からMSA表現に戻る情報フローを提供することでループを閉じ、全体的なEvoformerブロックがペア表現とMSA表現の間で情報を完全に混合し、構造モジュール内で構造生成の準備ができるようにする。

## End-to-end structure prediction
    structure module(構造予測head)はペア表現と, MSA表現のうち元の配列に対応するもののみを用いる。3次元骨格構造は, 各残基の独立な回転・並進移動により表される。この回転と並進は, 骨格部分の方向を優先する(...それを合わせることを目的に調整される?)ので, N-Calpha-Cに対する残基の方向は固定される。一方残基間の結合についての制約は無視されるので, 局所的な構造予測を優先するあまり残基の結合についての制約を無視することがよく観察される。Finetuning時にはこの構造違反に対する損失項を加えて制約を考慮させる。最終的には, 予測生成後にAmber force fieldのgradient descentという方法に沿って構造を微調整するときに制約が達成される。この微調整は構造予測精度を改善するわけではないが, 精度を損なわずに変な立体化学的違反を除去する(...立体化学的違反を除去する分精度が上がるが, 別の所で精度が落ちるので結局変わらないということ)?
The structure module (Fig. 3d) operates on a concrete 3D backbone structure using the pair representation and the original sequence row (single representation) of the MSA representation from the trunk. The 3D backbone structure is represented as Nres independent rotations and translations, each with respect to the global frame (residue gas) (Fig. 3e). These rotations and translations —representing the geometry of the N-Cα-C atoms— prioritize the orientation of the protein backbone so that the location of the side chain of each residue is highly constrained within that frame. Conversely, the peptide bond geometry is completely unconstrained and the network is observed to frequently violate the chain constraint during the application of the structure module as breaking this constraint enables the local refinement of all parts of the chain without solving complex loop closure problems. Satisfac- tion of the peptide bond geometry is encouraged during fine-tuning by a violation loss term. Exact enforcement of peptide bond geometry is only achieved in the post-prediction relaxation of the structure by gradient descent in the Amber32 force field. Empirically, this final relaxa- tion does not improve the accuracy of the model as measured by the global distance test (GDT)33 or lDDT-Cα34 but does remove distracting stereochemical violations without the loss of accuracy.
    構造モジュール（図3d）は、トランクからのMSA表現のペア表現と元の配列行（単一表現）を用いて、具体的な3Dバックボーン構造を操作する。3Dバックボーン構造は、Nres個の独立した回転と平行移動として表現され、それぞれがグローバルフレーム（残基ガス）を基準としている（図3e）。これらの回転と平行移動はN-Cα-C原子の形状を表し、タンパク質の背骨の向きを優先させるので、各残基の側鎖の位置はそのフレーム内で高度に拘束されます。逆に、ペプチド結合のジオメトリは全く制約されず、ネットワークは、複雑なループ閉鎖問題を解くことなく、鎖のすべての部分の局所的な精密化を可能にするため、構造モジュールの適用中に鎖の制約に頻繁に違反することが観察されます。ペプチド結合の形状を満足させることは、違反損失項によって微調整中に奨励される。ペプチド結合形状の厳密な強制は、Amber32力場の勾配降下による構造予測後の緩和でのみ達成される。経験的に、この最終的な緩和は、グローバル距離テスト（GDT）33やlDDT-Cα34によって測定されるモデルの精度を向上させないが、精度を損なうことなく、気が散るような立体化学的違反を除去する。

    structure moduleでは, まず配列表現を更新し, それをもとに構造を予測する。配列表現の更新にはinvariant point attention(IPA)というものを使う。これはMSAから計算したq, kによるattention weight, ペア表現のprojection weight, 前のモジュールで予測された座標と回転から計算されたattention weightを足し合わせ, それとMSA, 座標と回転のそれぞれから計算したvをかけ, さらにペア表現とのそれぞれの積の和?っぽいものを取り, それを足し合わせて新しい配列表現とする。その後, 配列表現から座標とquartanion(回転を4次元ベクトルで表現する手法。詳細：https://qiita.com/drken/items/0639cf34cce14e8d58a5)を予測, quartanionから回転行列を計算する。
    これは8層あるが, 重みを共有している。最終層以外の座標・回転予測は何をターゲットにして予測しているのかと思ったが重みを共有しているからいいのか。
The residue gas representation is updated iteratively in two stages (Fig. 3d). First, a geometry-aware attention operation that we term ‘invariant point attention’ (IPA) is used to update an Nres set of neural activations (single representation) without changing the 3D positions, then an equivariant update operation is performed on the residue gas using the updated activations. The IPA augments each of the usual attention queries, keys and values with 3D points that are produced in the local frame of each residue such that the final value is invariant to global rotations and translations (see Methods ‘IPA’ for details). The 3D queries and keys also impose a strong spatial/locality bias on the attention, which is well-suited to the iterative refinement of the protein structure. After each attention operation and element-wise transition block, the module computes an update to the rotation and translation of each backbone frame. The application of these updates within the local frame of each residue makes the overall attention and update block an equivariant operation on the residue gas.
    構造モジュール（図3d）は、トランクからのMSA表現のペア表現と元の配列行（単一表現）を用いて、具体的な3Dバックボーン構造を操作する。3Dバックボーン構造は、Nres個の独立した回転と平行移動として表現され、それぞれがグローバルフレーム（残基ガス）を基準としている（図3e）。これらの回転と平行移動はN-Cα-C原子の形状を表し、タンパク質の背骨の向きを優先させるので、各残基の側鎖の位置はそのフレーム内で高度に拘束されます。逆に、ペプチド結合のジオメトリは全く制約されず、ネットワークは、複雑なループ閉鎖問題を解くことなく、鎖のすべての部分の局所的な精密化を可能にするため、構造モジュールの適用中に鎖の制約に頻繁に違反することが観察されます。ペプチド結合の形状の満足は、残基ガス表現が2つの段階で反復的に更新される間に促されます（図3d）。まず、私たちが「不変点注意」（IPA）と呼ぶ形状を考慮した注意操作が、3次元位置を変えることなく神経活性（単一表現）のNresセットを更新するために使用され、次に、更新された活性を使用して残基ガスに対して等変量更新操作が実行される。IPAは、最終的な値がグローバルな回転と平行移動に対して不変であるように、各残渣のローカルフレームで生成される3D点で通常の注意クエリー、キー、値のそれぞれを補強する（詳細は方法「IPA」を参照）。3Dクエリーとキーはまた、注意に強い空間的／局所的バイアスを課し、これはタンパク質構造の反復的精密化に適している。各注意操作と要素ごとの遷移ブロックの後、モジュールは各骨格の回転と並進の更新を計算する。各残基のローカルフレーム内でこれらの更新を適用することで、全体的な注意と更新ブロックは残基gas.tuning上の等変量操作となります。


Predictions of side-chain χ angles as well as the final, per-residue accuracy of the structure (pLDDT) are computed with small per-residue networks on the final activations at the end of the network. The estimate of the TM-score (pTM) is obtained from a pairwise error prediction that is computed as a linear projection from the final pair representation. The final loss (which we term the frame-aligned point error (FAPE)(Fig.3f)) compares the predicted atom positions to the true positions under many different alignments. For each alignment, defined by aligning the predicted frame (Rk, tk) to the corresponding true frame, we com- pute the distance of all predicted atom positions xi from the true atom positions. The resulting Nframes × Natoms distances are penalized with a clamped L1 loss. This creates a strong bias for atoms to be correct relative to the local frame of each residue and hence correct with respect to its side-chain interactions, as well as providing the main source of chirality for AlphaFold (Supplementary Methods 1.9.3 and Supplementary Fig. 9).
    側鎖χ角の予測、および構造の最終的な残基あたりの精度(pLDDT)は、小さな残基あたりネットワークを用いて、ネットワークの末端での最終的な活性化で計算される。TMスコアの推定値（pTM）は、最終的なペア表現からの線形射影として計算されるペアごとの誤差予測から得られる。finalalloss(これはframe-alignedpointerror(FAPE)(図3f)を決定する)は、多くの異なるアライメントの下で予測された原子位置と真の位置を比較する。予測フレーム(Rk, tk)を対応する真のフレームにアライメントすることで定義される各アライメントについて、全ての予測原子位置xiの真の原子位置からの距離を計算する。結果として得られるNframes×Natomsの距離はクランプL1ロスでペナルティを受ける。これにより、各残基のローカルフレームに対して原子が正しくなるように強いバイアスがかかり、その結果、側鎖相互作用に関して正しくなると同時に、AlphaFoldのキラリティの主要なソースとなる（補足方法1.9.3および補足図9）。

## Training with labelled and unlabelled data
一度PDBだけでモデルを学習させた後, Uniclust30という35万の立体構造のないタンパク質のアミノ酸配列のデータセットについて構造を予測する。その中でconfidenceの高かったものを抽出し, それとPDBを混ぜ, もう一度同じモデルをscratchから学習させる。(アミノ酸配列の?)croppingなどaugmentationをかけることで, 同じ構造を簡単に予測できないようにする。これによって精度が大幅に高まった。

    AlphaFoldは, PDBからの教師あり学習のみからでも高い精度を
The AlphaFold architecture is able to train to high accuracy using only supervised learning on PDB data, but we are able to enhance accuracy (Fig. 4a) using an approach similar to noisy student self-distillation35. In this procedure, we use a trained network to predict the structure of around 350,000 diverse sequences from Uniclust3036 and make a new dataset of predicted structures filtered to a high-confidence subset. We then train the same architecture again from scratch using a mixture of PDB data and this new dataset of predicted structures as the training data, in which the various training data augmentations such as crop- ping and MSA subsampling make it challenging for the network to recapitulate the previously predicted structures. This self-distillation procedure makes effective use of the unlabelled sequence data and considerably improves the accuracy of the resulting network.
    AlphaFoldアーキテクチャは、PDBデータに対して教師あり学習のみを用いて高精度に学習することができるが、ノイズの多い学生の自己蒸留35に似たアプローチを用いて精度を高めることができる（図4a）。この手順では、Uniclust3036から約35万個の多様な配列の構造を予測するために訓練されたネットワークを使用し、高信頼性のサブセットにフィルタリングされた予測構造の新しいデータセットを作成する。その際、クロッピングやMSAサブサンプリングなどの様々な学習データの増強により、ネットワークが以前に予測した構造を再現することが困難になる。この自己蒸留法は、ラベル付けされていない配列データを有効に利用し、得られたネットワークの精度を大幅に向上させる。

Additionally, we randomly mask out or mutate individual residues within the MSA and have a Bidirectional Encoder Representations from Transformers (BERT)-style37 objective to predict the masked elements of the MSA sequences. This objective encourages the network to learn to interpret phylogenetic and covariation relationships without hardcoding a particular correlation statistic into the features. The BERT objective is trained jointly with the normal PDB structure loss on the same training examples and is not pre-trained, in contrast to recent independent work38.
    さらに、MSA内の個々の残基をランダムにマスクアウトまたは変異させ、MSA配列のマスクされた要素を予測するBERT(Bidirectional Encoder Representations from Transformers)スタイルの目的37を持つ。この目的は、特定の相関統計量を特徴量にハードコーディングすることなく、系統関係と共分散関係を解釈する学習をネットワークに促す。BERT目的は、同じ訓練例で通常のPDB構造損失と共同で訓練され、最近の独立した研究38とは対照的に、事前訓練されていない。

## Interpreting the neural network
Alphafold2がどのように構造を予測するのか調べるため, Evoformerをブロックごとに学習させた(...どういうこと?)
最初の数層でおおまかな構造を作り, それを少しずつ改善している感じだった。
SARS-CoV-2のなんとかとかいうタンパク質など, 難しい構造については何度も構造を配置し直していた。


To understand how AlphaFold predicts protein structure, we trained a separate structure module for each of the 48 Evoformer blocks in the network while keeping all parameters of the main network fro- zen (Supplementary Methods 1.14). Including our recycling stages, this provides a trajectory of 192 intermediate structures—one per full Evoformer block—in which each intermediate represents the belief of the network of the most likely structure at that block. The resulting trajectories are surprisingly smooth after the first few blocks, show- ing that AlphaFold makes constant incremental improvements to the structure until it can no longer improve (see Fig. 4b for a trajectory of accuracy). These trajectories also illustrate the role of network depth. For very challenging proteins such as ORF8 of SARS-CoV-2 (T1064), the network searches and rearranges secondary structure elements for many layers before settling on a good structure. For other proteins such as LmrP (T1024), the network finds the final structure within the first few layers. Structure trajectories of CASP14 targets T1024, T1044, T1064 and T1091 that demonstrate a clear iterative building process for a range of protein sizes and difficulties are shown in Supplementary Videos 1–4. In Supplementary Methods 1.16 and Supplementary Figs. 12, 13, we interpret the attention maps produced by AlphaFold layers.
    AlphaFoldがどのようにタンパク質の構造を予測するのかを理解するために、主要なネットワークのパラメータはすべてそのままに、ネットワーク内の48のEvoformerブロックごとに個別の構造モジュールをトレーニングした（補足方法1.14）。リサイクル段階を含めると、192の中間構造の軌跡が得られ、各中間構造はそのブロックで最も可能性の高い構造のネットワークの信念を表している。得られた軌跡は、最初の数ブロック以降は驚くほど滑らかであり、AlphaFoldがもはや改善できなくなるまで、構造に対して絶え間ない漸進的な改善を行うことを示している（精度の軌跡については図4bを参照）。これらの軌跡は、ネットワークの深さの役割も示している。SARS-CoV-2のORF8（T1064）のような非常に困難なタンパク質では、ネットワークは、良好な構造に落ち着く前に、何層にもわたって二次構造要素を検索し、再配置します。LmrP（T1024）のような他のタンパク質では、ネットワークは最初の数層で最終構造を見つける。CASP14のターゲットであるT1024, T1044, T1064, T1091の構造軌跡は、様々なタンパク質のサイズと難易度に対する明確な反復構築プロセスを示しており、補足ビデオ1-4に示されている。補足方法1.16と補足図12、13では、AlphaFold層によって生成されたアテンションマップを解釈しています。

Figure 4a contains detailed ablations of the components of AlphaFold that demonstrate that a variety of different mechanisms contribute to AlphaFold accuracy. Detailed descriptions of each ablation model, their training details, extended discussion of ablation results and the effect of MSA depth on each ablation are provided in Supplementary Methods 1.13 and Supplementary Fig. 10.
    図4aには、AlphaFoldの構成要素の詳細なアブレーションが含まれており、様々な異なるメカニズムがAlphaFoldの精度に寄与していることが示されている。各アブレーション・モデルの詳細な説明、トレーニングの詳細、アブレーション結果の拡大考察、各アブレーションに対するMSA深度の影響は、補足方法1.13と補足図10に記載されている。

## MSA depth and cross-chain contacts
AlphaFoldの欠点としては, まずはMSAのalignment depth(よくわからないけど数ということ?)が30以下だと精度が落ちる
他のタンパク質の間を取り持つもので, その構造が他のタンパク質によって決まるようなものについても精度が低い。
逆にホモマーならば, 同じタンパク質同士が複雑に結合していても予測できる。

Although AlphaFold has a high accuracy across the vast majority of deposited PDB structures, we note that there are still factors that affect accuracy or limit the applicability of the model. The model uses MSAs and the accuracy decreases substantially when the median alignment depth is less than around 30 sequences (see Fig. 5a for details). We observe a threshold effect where improvements in MSA depth over around 100 sequences lead to small gains. We hypothesize that the MSA information is needed to coarsely find the correct structure within the early stages of the network, but refinement of that prediction into a high-accuracy model does not depend crucially on the MSA information. The other substantial limitation that we have observed is that AlphaFold is much weaker for proteins that have few intra-chain or homotypic con- tacts compared to the number of heterotypic contacts (further details are provided in a companion paper39). This typically occurs for bridging domains within larger complexes in which the shape of the protein is created almost entirely by interactions with other chains in the complex. Conversely, AlphaFold is often able to give high-accuracy predictions for homomers, even when the chains are substantially intertwined (Fig. 5b). We expect that the ideas of AlphaFold are readily applicable to predicting full hetero-complexes in a future system and that this will remove the dif- ficulty with protein chains that have a large number of hetero-contacts.
    AlphaFoldは、寄託されたPDB構造の大部分において高い精度を示していますが、精度に影響を与えたり、モデルの適用性を制限したりする要因がまだあることに注意してください。このモデルはMSAを使用しており、アライメントの深さの中央値が約30配列未満になると精度が大幅に低下します（詳細は図5aを参照）。また、MSAの深さが100配列以上になると、閾値効果が現れ、精度が向上する。我々は、MSA情報はネットワークの初期段階で正しい構造を粗く見つけるために必要であるが、その予測を高精度モデルに洗練することは、MSA情報に大きく依存しないという仮説を立てた。我々が観察したもう一つの重要な限界は、AlphaFoldは、ヘテロ型コンタクトの数に比べて、鎖内またはホモ型コンタクトが少ないタンパク質に対しては、かなり弱いということである（さらなる詳細は、関連論文に記載されている39）。これは、タンパク質の形状がほとんど複合体中の他の鎖との相互作用によって作られるような、大きな複合体中の橋渡しドメインで典型的に起こる。逆にAlphaFoldは、鎖がかなり絡み合っていても、ホモマーに対して高精度の予測を行うことができる場合が多い（図5b）。AlphaFoldのアイデアは、将来のシステムにおいて、完全なヘテロ複合体の予測に容易に適用でき、多数のヘテロ接点を持つタンパク質鎖の困難さを取り除くことができると期待している。

## Related work
The prediction of protein structures has had a long and varied develop- ment, which is extensively covered in a number of reviews14,40–43. Despite the long history of applying neural networks to structure prediction14,42,43, they have only recently come to improve structure prediction10,11,44,45. These approaches effectively leverage the rapid improvement in com- puter vision systems46 by treating the problem of protein structure prediction as converting an ‘image’ of evolutionary couplings22–24 to an ‘image’ of the protein distance matrix and then integrating the distance predictions into a heuristic system that produces the final 3D coordinate prediction. A few recent studies have been developed to predict the 3D coordinates directly47–50, but the accuracy of these approaches does not match traditional, hand-crafted structure prediction pipelines51. In paral- lel, the success of attention-based networks for language processing52 and, more recently, computer vision31,53 has inspired the exploration of attention-based methods for interpreting protein sequences54–56.
    タンパク質の構造予測には長い歴史と様々な発展があり、多くの総説14,40-43がある。ニューラルネットワークを構造予測に応用した歴史は長いにもかかわらず14,42,43、構造予測を改善するようになったのはごく最近のことである10,11,44,45。これらのアプローチは、進化的カップリング22-24の「画像」をタンパク質の距離行列の「画像」に変換し、距離予測をヒューリスティックシステムに統合して最終的な3次元座標予測を生成するものとしてタンパク質構造予測の問題を扱うことにより、コンピュータビジョンシステム46の急速な改善を効果的に利用している。最近、3次元座標を直接予測する研究がいくつか開発された47-50が、これらのアプローチの精度は、従来の手作りの構造予測パイプライン51には及ばない。パラレルでは、言語処理52や、最近ではコンピュータビジョン31,53における注意力ベースのネットワークの成功が、タンパク質配列の解釈のための注意力ベースの手法の探求を促している54-56。


## Discussion
　AlphaFoldは, 手動でいろんな相互作用を設定することなくそれらを学習しているのがいいよね(例えば, 水素結合について設定しなくとも水素結合を生成する)
The methodology that we have taken in designing AlphaFold is a combi- nation of the bioinformatics and physical approaches: we use a physical and geometric inductive bias to build components that learn from PDB data with minimal imposition of handcrafted features (for example, AlphaFold builds hydrogen bonds effectively without a hydrogen bond score function). This results in a network that learns far more efficiently from the limited data in the PDB but is able to cope with the complexity and variety of structural data.
    例えば、AlphaFoldは水素結合スコア関数を使わずに水素結合を効率的に構築する。この結果、PDBの限られたデータからはるかに効率的に学習するネットワークが得られるが、構造データの複雑さと多様性に対応することができる。

In particular, AlphaFold is able to handle missing the physical context and produce accurate models in challenging cases such as intertwined homomers or proteins that only fold in the presence of an unknown haem group. The ability to handle underspecified structural conditions is essential to learning from PDB structures as the PDB represents the full range of conditions in which structures have been solved. In gen- eral, AlphaFold is trained to produce the protein structure most likely to appear as part of a PDB structure. For example, in cases in which a particular stochiometry, ligand or ion is predictable from the sequence alone, AlphaFold is likely to produce a structure that respects those constraints implicitly.
    特にAlphaFoldは、物理的な文脈が欠落していても扱うことができ、絡み合ったホモマーや未知のヘム基の存在下でのみ折り畳まれるタンパク質などの困難なケースでも正確なモデルを生成することができる。PDBは構造が解かれたあらゆる条件を表しているため、指定されていない構造条件を扱う能力は、PDB構造から学習するのに不可欠である。通常、AlphaFoldはPDB構造の一部として現れる可能性が最も高いタンパク質構造を生成するように学習される。例えば、特定のストキオメトリー、リガンド、イオンが配列のみから予測可能な場合、AlphaFoldは暗黙のうちにそれらの制約を尊重した構造を生成する可能性が高い。

AlphaFoldは1タンパク質あたり数分~数時間なので, すぐ終わる


AlphaFold has already demonstrated its utility to the experimental community, both for molecular replacement57 and for interpreting cryogenic electron microscopy maps58. Moreover, because AlphaFold outputs protein coordinates directly, AlphaFold produces predictions in graphics processing unit (GPU) minutes to GPU hours depending on the length of the protein sequence (for example, around one GPU min- ute per model for 384 residues; see Methods for details). This opens up the exciting possibility of predicting structures at the proteome-scale and beyond—in a companion paper39, we demonstrate the application of AlphaFold to the entire human proteome39.
    AlphaFoldは、分子置換57と低温電子顕微鏡マップ58の解釈の両方において、すでに実験コミュニティでその有用性が実証されている。さらに、AlphaFoldはタンパク質の座標を直接出力するため、タンパク質配列の長さにもよるが、GPU（グラフィック・プロセッシング・ユニット）分～GPU時間で予測結果を得ることができる（例えば、384残基の場合、1モデルあたり約1GPU分。） これは、プロテオームスケールでの構造予測や、それ以上の構造予測の可能性を開くものである。

The explosion in available genomic sequencing techniques and data has revolutionized bioinformatics but the intrinsic challenge of experi- mental structure determination has prevented a similar expansion in our structural knowledge. By developing an accurate protein structure prediction algorithm, coupled with existing large and well-curated structure and sequence databases assembled by the experimental community, we hope to accelerate the advancement of structural bioinformatics that can keep pace with the genomics revolution. We hope that AlphaFold—and computational approaches that apply its techniques for other biophysical problems—will become essential tools of modern biology.
    利用可能なゲノム配列決定技術とデータの爆発的な増加は、バイオインフォマティクスに革命をもたらしたが、実験的構造決定という本質的な課題が、構造知識の同様の拡大を妨げてきた。私たちは、正確なタンパク質構造予測アルゴリズムを開発し、実験コミュニティによって構築された既存の大規模かつ十分にキュレーションされた構造・配列データベースと組み合わせることで、ゲノム革命と歩調を合わせることができる構造バイオインフォマティクスの進歩を加速させたいと考えています。我々は、AlphaFoldと、その技術を他の生物物理学的問題に応用する計算機的アプローチが、現代生物学の必須ツールとなることを願っている。
