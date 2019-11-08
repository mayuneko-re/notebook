<!---
  コメント
-->

理論は教科書に書いてある．
自分でやってみればもっと分かるはず．
油層工学に関する面白そうな話を Jupyter notebooks (Google Colaboratory) で実装してみる．


## Buckley-Leverett Solution

- [水油置換問題](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Basic.ipynb)
- [水油置換問題（感度解析）](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Basic_Sensitivity.ipynb)
- [水油置換問題（重力の影響込み）](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Gravity.ipynb)

## フラッシュ計算

- [PTフラッシュ計算](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/PT_Flash_Calculation.ipynb)
- [Ternary Phase Diagram](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Ternary_Phase_Diagram.ipynb)



## 参考文献

*   Stalkup, F. I. (1983). *Miscible displacement*.
*   Lake, L. W. (1989). *Enhanced oil recovery*. Old Tappan, NJ: Prentice Hall Inc.
*   Orr, F. M. (2007). *Theory of gas injection processes*. Copenhagen: Tie-Line Publications.
*   [PetroWiki](https://petrowiki.org/PetroWiki)


<ul>
  {% for post in site.posts %}
    <li>
      <a href="{{ post.url }}">{{ post.title }}</a>
    </li>
  {% endfor %}
</ul>

<ul class="posts">
  {% for post in site.posts %}
    <li>
      <span>{{ post.date | date_to_string }}</span> » <a href="{{ post.url | relative_url }}" title="{{ post.title }}">{{ post.title }}</a>
    </li>
  {% endfor %}
</ul>
