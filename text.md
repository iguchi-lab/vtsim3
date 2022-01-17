# 温熱計算プログラムvtsim3の使い方

---
##1.はじめに

　vtsimは、建築環境工学・建築設備を学ぶ大学生、研究者向けに開発された温熱計算プログラムで下記の特徴があります。

*   ブラウザ上（Google Clab）で動作する。環境設定が不要。
*   自由度の高いプログラム言語pythonで記述。
*   速度が求められる計算部分はc++で記述。
*   ユーザーは、プログラム言語pythonに加え、pandas、numpyなどの知識が必要。
*   熱・換気回路網による計算をベースとし、節点（ノード）と回路網（ネットワーク）、各種条件の設定で動作。


##2.インストール方法

　Google Colabで下記のコードを実行して、インストールする。

```
!pip install git+https://github.com/iguchi-lab/vtsim3
from vtsim import vtsim as vt

```
　あわせて下記も実行しておくとよい。

```
!pip install japanize-matplotlib
import matplotlib.pyplot as plt
import japanize_matplotlib

import numpy as np
import pandas as pd
```
