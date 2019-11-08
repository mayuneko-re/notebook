理論は教科書に書いてある．
自分でやってみればもっと分かるはず．
油層工学に関する面白そうな話を Jupyter notebooks (Google Colaboratory) で実装してみる．


## Buckley-Leverett Solution

- [水油置換問題 :octocat:](/colab/Buckley_Leverett_Basic_Sensitivity.ipynb)
[:page_facing_up:](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Basic.ipynb)
- [水油置換問題（感度解析） :octocat:](/colab/Buckley_Leverett_Basic_Sensitivity.ipynb)
[:page_facing_up:](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Basic_Sensitivity.ipynb)
- [水油置換問題（重力の影響込み） :octocat:](/colab/Buckley_Leverett_Gravity.ipynb)
[:page_facing_up:](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Buckley_Leverett_Gravity.ipynb)

## フラッシュ計算

- [PTフラッシュ計算 :octocat:](/colab/PT_Flash_Calculation.ipynb)
[:page_facing_up:](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/PT_Flash_Calculation.ipynb)
- [Ternary Phase Diagram :octocat:](/colab/Ternary_Phase_Diagram.ipynb)
[:page_facing_up:](https://nbviewer.jupyter.org/github/mayuneko-re/notebook/blob/master/colab/Ternary_Phase_Diagram.ipynb)



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
