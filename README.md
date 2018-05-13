# ModelRocketSimulator

OverView \\
モデルロケットの数値シミュレーションです．特別なことは行っていません．回転の運動方程式で角速度を出した後，その角速度にそって次のクォータニオンを求める．そのクォータニオンを使用して機軸を変更して並進運動方程式を解いています．

## Demo
MRS_demo.mv4をご覧ください．

## Requirement
numpy \\
matplotlib \\
pandas \\
scipy \\
Anacondaを入れれば一発のはずです

## Usage
・spec.csvで機体スペックの調整 \\
・G75.txtは G75というモデロケのエンジンの推力履歴です.　他の推力データは http://www.thrustcurve.org からダウンロードしましょう．\\
・MRS.pyを実行すればz-tグラフ，x,y,zの3次元グラフが見れます．

## Install
