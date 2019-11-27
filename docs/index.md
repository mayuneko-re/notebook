理論は教科書に書いてある．
それじゃあやってみよう，ということで，油層工学に関する面白そうな話を実装してみるサイト．
興味は，多相流動とか，EORとか，もっと言えばCO2-EORとか．

<ul class="posts">
  {% for post in site.posts %}
    <li>
      <small>{{ post.date | date: "%-d %B %Y" }}</small><br>
      <a href="{{ post.url | relative_url }}" title="{{ post.title }}">{{ post.title }}</a>
    </li>
  {% endfor %}
</ul>


