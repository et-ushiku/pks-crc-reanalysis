"""CN-CRC（Yu 2017、中国、PRJEB10878）128検体を ENA から落とす。

AT-CRC 版（at_download.py）からの変更点
  - URL は cn_pick.tsv の4列目に入れておいたので、ENA の表を別に読まない
  - 1検体のファイルは **2個**（AT-CRC は3個）。接続9本を埋めるには
    落としかけの検体を多めに持つ必要があるので MAX_INFLIGHT を 6 にした
  - 検体名は cancer_ERRxxxxxxx / n_ERRxxxxxxx（cn_pick.tsv の1列目）をそのまま使う

速度について（AT-CRC での実測、2026-11-16）
  ENA は接続1本あたり約0.5〜0.7 MB/s に絞っている。回線自体は余っていて、
  本数を増やすと1本の速度は落ちずに合計だけ増えた。
  ※ ENA は公的資金で運営される共有資源なので、本数は上げすぎない（9本に留める）。

  610.4 GB を 4.75 MB/s で約 37 時間の見込み。

流し込み方式: 落とす → 当てる(cn_map.sh) → 消す を並行させ、ディスクを使いすぎない。
  - 1検体ぶんのファイルが揃ったら .complete を置く（cn_map.sh がこの印を待つ）
  - cn_map.sh が終わると .mapped を置き、fastq を消す
  - 未処理の検体が溜まったら、落とすのを一時停止する

使い方: python3 cn_download.py
"""
import os
import queue
import subprocess
import threading
import time

META = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(os.path.dirname(META), 'data_cn')
PICK = os.path.join(META, 'cn_pick.tsv')

WORKERS = 9       # 同時に張る接続の数（常にこれだけ埋める）
MAX_INFLIGHT = 6  # 同時に落としかけの検体数（1検体2ファイルなので多めに持つ）
AHEAD = 3         # 当て終わっていない検体がこれ以上溜まったら落とすのを止める

_lock = threading.Lock()


def log(msg):
    with _lock:
        print(f'{time.strftime("%m/%d %H:%M:%S")} {msg}', flush=True)


def valid_gzip(path):
    return subprocess.run(['gzip', '-t', path], capture_output=True).returncode == 0


def fetch(url_path, dest):
    """通信が切れても -C - で続きから取る。中身が壊れていたときだけ消して取り直す。"""
    if os.path.exists(dest) and valid_gzip(dest):
        return True
    for attempt in range(1, 21):
        proto = 'https' if attempt % 2 else 'ftp'
        r = subprocess.run(['curl', '-sS', '-f', '-C', '-', '--connect-timeout', '60',
                            '--speed-limit', '20000', '--speed-time', '300',
                            '-o', dest, f'{proto}://{url_path}'],
                           capture_output=True, text=True)
        if r.returncode == 0:
            if valid_gzip(dest):
                return True
            log(f'  壊れたファイル、消して取り直し {os.path.basename(dest)}')
            os.remove(dest)
        else:
            size = os.path.getsize(dest) / 1e6 if os.path.exists(dest) else 0
            log(f'  中断 ({attempt}回目, {proto}, {size:.0f}MBまで) '
                f'{os.path.basename(dest)}: {r.stderr.strip()[:100]}')
        time.sleep(min(30 * attempt, 300))
    return False


def pending():
    """落とし終わったが、まだ当て終わっていない検体の数"""
    if not os.path.isdir(BASE):
        return 0
    return sum(1 for d in os.listdir(BASE)
               if os.path.exists(f'{BASE}/{d}/.complete')
               and not os.path.exists(f'{BASE}/{d}/.mapped'))


class Sample:
    def __init__(self, idx, name, d, files):
        self.idx, self.name, self.d, self.files = idx, name, d, files
        self.left = len(files)
        self.ok = True
        self.t0 = time.time()


def main():
    picks = []
    with open(PICK, encoding='utf-8') as f:
        next(f)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 4:
                picks.append((p[0], [u for u in p[3].split(';') if u]))

    os.makedirs(BASE, exist_ok=True)

    # まだ落とし終わっていない検体だけを、元の順（がん/対照を比率どおりに混ぜた順）で拾う
    todo = []
    for i, (name, files) in enumerate(picks, 1):
        d = f'{BASE}/{name}'
        if os.path.exists(f'{d}/.mapped') or os.path.exists(f'{d}/.complete'):
            continue
        if not files:
            log(f'[{i}/{len(picks)}] URL がない {name}')
            continue
        todo.append(Sample(i, name, d, files))
    log(f'残り {len(todo)} 検体 / 全 {len(picks)} 検体。接続 {WORKERS} 本を常に埋める')

    q = queue.Queue()
    stats = {'bytes': 0, 'done': 0}
    t_start = time.time()

    def worker():
        while True:
            item = q.get()
            if item is None:
                q.task_done()
                return
            smp, url = item
            dest = f'{smp.d}/{url.rsplit("/", 1)[-1]}'
            ok = fetch(url, dest)
            with _lock:
                if not ok:
                    smp.ok = False
                else:
                    stats['bytes'] += os.path.getsize(dest)
                smp.left -= 1
                finished = smp.left == 0
            if finished:
                if smp.ok:
                    mb = sum(os.path.getsize(f'{smp.d}/{x}')
                             for x in os.listdir(smp.d) if x.endswith('.gz')) / 1e6
                    open(f'{smp.d}/.complete', 'w').close()
                    with _lock:
                        stats['done'] += 1
                        n = stats['done']
                        rate = stats['bytes'] / 1e6 / (time.time() - t_start)
                    log(f'[{smp.idx}/{len(picks)}] 完了 {smp.name} {mb:.0f}MB '
                        f'{(time.time()-smp.t0)/60:.1f}分  '
                        f'累計 {n} 検体 平均 {rate:.2f} MB/s')
                else:
                    open(f'{smp.d}/.failed', 'w').close()
                    log(f'[{smp.idx}/{len(picks)}] 失敗 {smp.name}（再実行すれば続きから）')
            q.task_done()

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(WORKERS)]
    for t in threads:
        t.start()

    # 落としかけの検体を MAX_INFLIGHT 個までに保ちながら、順に投入していく
    it = iter(todo)
    inflight = []
    exhausted = False
    while not exhausted or inflight:
        inflight = [s for s in inflight if s.left > 0]
        while len(inflight) < MAX_INFLIGHT and not exhausted:
            if pending() >= AHEAD:
                break
            try:
                smp = next(it)
            except StopIteration:
                exhausted = True
                break
            os.makedirs(smp.d, exist_ok=True)
            inflight.append(smp)
            log(f'[{smp.idx}/{len(picks)}] 開始 {smp.name} ({len(smp.files)}ファイル)')
            for u in smp.files:
                q.put((smp, u))
        time.sleep(20)

    q.join()
    for _ in threads:
        q.put(None)

    el = time.time() - t_start
    log(f'全検体のダウンロード処理が終了。{stats["bytes"]/1e9:.1f}GB を {el/3600:.1f}時間 '
        f'= 平均 {stats["bytes"]/1e6/el:.2f} MB/s')
    open(f'{BASE}/.download_finished', 'w').close()


if __name__ == '__main__':
    main()
