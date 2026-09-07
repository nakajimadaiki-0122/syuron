"""
plotstyle.py -- 図の共通設定
============================

matplotlib の既定フォント（DejaVu Sans）には日本語が無く、
図中の日本語が豆腐（□）になる。Windows にある日本語フォントへ切り替える。

    import plotstyle; plotstyle.use_jp()
"""
import matplotlib
import matplotlib.font_manager as fm

CANDIDATES = ["Yu Gothic", "Meiryo", "MS Gothic", "MS UI Gothic",
              "Noto Sans CJK JP", "Hiragino Sans", "IPAexGothic"]


def use_jp():
    """使える日本語フォントを探して設定する。見つかった名前を返す。"""
    have = {f.name for f in fm.fontManager.ttflist}
    for name in CANDIDATES:
        if name in have:
            matplotlib.rcParams["font.family"] = name
            matplotlib.rcParams["axes.unicode_minus"] = False
            return name
    return None
