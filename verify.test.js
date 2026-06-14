// PSU（電源）と RD（反応拡散）の核心式を、実装からそのまま写して検証する。推測なし。
let pass=0,fail=0,log=[];
function ok(n,c,x){ if(c){pass++;log.push('  ok  '+n);} else {fail++;log.push('  FAIL '+n+(x?'  '+x:''));} }
const close=(a,b,e=1e-6)=>Math.abs(a-b)<=e*(1+Math.abs(b));

// ===== PSU: operatingPoint を写経 =====
function operatingPoint(p){
  const {topology,Vin,Vref:Vout,Iout,Rload:R,L,rL,Cout,ESR,fsw}=p;
  if(topology==="boost"){
    let Dp=Vin/Vout; Dp=Math.min(Math.max(Dp,1e-3),0.999); const D=1-Dp;
    const Pout=Vout*Iout,Iin=Pout/Vin,ILavg=Iin;
    const dIL=(Vin*D)/(L*fsw);
    const f_lc=Dp/(2*Math.PI*Math.sqrt(L*Cout));
    const f_rhp=(Dp*Dp*R)/(2*Math.PI*L);
    const mode=ILavg<dIL/2?"DCM":"CCM";
    return {D,Pout,Iin,dIL,f_lc,f_rhp,mode};
  }
  const D=Vout/Vin,Pout=Vout*Iout,Iin=Pout/Vin;
  const dIL=((Vin-Vout)*D)/(L*fsw);
  const f_lc=1/(2*Math.PI*Math.sqrt(L*Cout));
  const mode=Iout<dIL/2?"DCM":"CCM";
  return {D,Pout,Iin,dIL,f_lc,f_rhp:Infinity,mode};
}
const base=(o)=>Object.assign({topology:"buck",Vin:12,Vref:6,Iout:2,Rload:5,L:47e-6,rL:0.03,Cout:220e-6,ESR:0.02,fsw:3e5},o);

// Buck
let b=operatingPoint(base({topology:"buck",Vin:12,Vref:6}));
ok('Buck Vin=12,Vout=6 → D=0.5', close(b.D,0.5), 'D='+b.D);
ok('Buck f_rhp=∞（RHP零点なし）', b.f_rhp===Infinity);
let b2=operatingPoint(base({topology:"buck",Vin:12,Vref:9}));
ok('Buck Vout↑ で D↑', b2.D>b.D, b2.D+'>'+b.D);
let bL1=operatingPoint(base({topology:"buck",L:47e-6})), bL2=operatingPoint(base({topology:"buck",L:200e-6}));
ok('Buck L↑ でインダクタリップル↓', bL2.dIL<bL1.dIL);
ok('Buck Pout=Vout*Iout', close(b.Pout,6*2));

// Boost
let g=operatingPoint(base({topology:"boost",Vin:5,Vref:10,Rload:5,L:150e-6}));
ok('Boost Vin=5,Vout=10 → D=0.5', close(g.D,0.5), 'D='+g.D);
ok('Boost f_rhp 有限', isFinite(g.f_rhp));
// f_rhp = Dp^2 R/(2πL), Dp=0.5, R=5, L=150e-6 → 0.25*5/(2π*150e-6)
let expect=(0.5*0.5*5)/(2*Math.PI*150e-6);
ok('Boost f_rhp 手計算一致', close(g.f_rhp,expect,1e-9), g.f_rhp+' vs '+expect);
// 単調性: R↓→f_rhp↓, L↑→f_rhp↓, D↑(高昇圧比)→f_rhp↓
let gR=operatingPoint(base({topology:"boost",Vin:5,Vref:10,Rload:2,L:150e-6}));
ok('Boost R↓ で f_RHP↓', gR.f_rhp<g.f_rhp, gR.f_rhp+'<'+g.f_rhp);
let gL=operatingPoint(base({topology:"boost",Vin:5,Vref:10,Rload:5,L:400e-6}));
ok('Boost L↑ で f_RHP↓', gL.f_rhp<g.f_rhp);
let gD=operatingPoint(base({topology:"boost",Vin:3,Vref:15,Rload:5,L:150e-6})); // 高昇圧比 D大
ok('Boost 高昇圧比(D大) で f_RHP↓', gD.D>g.D && gD.f_rhp<g.f_rhp, 'D='+gD.D+' frhp='+gD.f_rhp);

// Bode: 20log10
const mag=(x)=>20*Math.log10(x);
ok('Bode 20log10(1)=0dB', close(mag(1),0));
ok('Bode 20log10(10)=20dB', close(mag(10),20));
ok('Bode 20log10(0.5)≈-6.02dB', close(mag(0.5),-6.0206,1e-3));

// ===== RD: 9点ステンシル＋オイラー1ステップを写経 =====
const wC=-1, wO=0.20, wD=0.05;
ok('RD ステンシル重み合計=0（定数場で∇²=0）', close(wC+4*wO+4*wD,0));
// 一様 U=1,V=0 はオイラー1ステップで不変
function rdStep(u,v,lapU,lapV,F,k,Du,Dv,dt){
  const uvv=u*v*v;
  let un=u+(Du*lapU - uvv + F*(1-u))*dt;
  let vn=v+(Dv*lapV + uvv - (F+k)*v)*dt;
  return [Math.min(Math.max(un,0),1), Math.min(Math.max(vn,0),1)];
}
// 一様場: lapU=lapV=0（重み合計0なので）
let [un,vn]=rdStep(1,0, 0,0, 0.0367,0.0649,0.16,0.08,1.0);
ok('RD 一様 U=1,V=0 は1ステップで不変', close(un,1)&&close(vn,0), 'U='+un+' V='+vn);
// 中央にVの種があると、隣のマスに正のlapVが効いてVが増える（拡散）
// 隣マス: 中心V=1, 自分V=0 → lapV = wO*1 + ... > 0
let lapV_neighbor = wO*1; // 直交方向に種が1つ隣接（簡略）
let [, vneigh]=rdStep(1,0, 0,lapV_neighbor, 0.0367,0.0649,0.16,0.08,1.0);
ok('RD 種の隣でVが拡散して増える(>0)', vneigh>0, 'V='+vneigh);

console.log(log.join('\n'));
console.log('\n==== '+pass+' passed, '+fail+' failed ====');
process.exit(fail?1:0);
