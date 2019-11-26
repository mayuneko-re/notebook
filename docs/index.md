理論は教科書に書いてある．
それじゃあやってみよう，ということで，油層工学に関する面白そうな話を実装してみるサイト．
興味は，多相流動とか，EORとか，もっと言えばCO2-EORとか．

## 油層ノート

- [Buckley-Leverettの解法で水油置換問題を解く](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Basic.ipynb)
- [Buckley-Leverettの解法で水油置換問題を解く（その２：粘性と総体浸透率の感度解析）](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Basic_Sensitivity.ipynb)
- [Buckley-Leverettの解法で水油置換問題を解く（その３：重力の影響込み）](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Gravity.ipynb)
- [Koval法でミシブル置換をモデルする](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Koval_method_for_miscible_displacement.ipynb)
- [Peng-Robinsonの状態方程式でPTフラッシュ計算をする](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/PT_Flash_Calculation.ipynb)
- [Peng-Robinsonの状態方程式でTernary Phase Diagramを描いてみる](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Ternary_Phase_Diagram.ipynb)
- [マルチセルモデルでCO2圧入に伴う油層流体の挙動を考える](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Multi_Cell_Model_for_CO2_core_flooding.ipynb)


## メモ等の置き場

<ul class="posts">
  {% for post in site.posts %}
    <li>
      <a href="{{ post.url | relative_url }}" title="{{ post.title }}">{{ post.title }}</a>
    </li>
  {% endfor %}
</ul>

## 参考文献

*   Stalkup, F. I. (1983). *Miscible displacement*.
*   Lake, L. W. (1989). *Enhanced oil recovery*. Old Tappan, NJ: Prentice Hall Inc.
*   Orr, F. M. (2007). *Theory of gas injection processes*. Copenhagen: Tie-Line Publications.
*   [PetroWiki](https://petrowiki.org/PetroWiki)

