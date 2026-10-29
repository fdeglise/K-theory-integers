
const output={};
for(const n of [8,9]) {
 const bound=n===8?19:37, byOrder={};
 for(let d=2;d<=bound;d++) {
  const h=Math.floor(d/2), weights=Array.from({length:h+1},(_,a)=>Array.from({length:h},(_,j)=>{let r=((j+1)*a)%d;return Math.min(r,d-r)}));
  const maxTail=Array.from({length:h+1},()=>Array(h).fill(0));
  for(let a=h;a>=0;a--) for(let j=0;j<h;j++)maxTail[a][j]=Math.max(weights[a][j],a<h?maxTail[a+1][j]:0);
  const gcd=(a,b)=>b?gcd(b,a%b):a;
  const units=[];for(let u=2;u<=h;u++)if(gcd(u,d)===1)units.push(u);
  const candidates=[],a=[],sums=Array(h).fill(0);
  const lexless=(a,b)=>{for(let i=0;i<a.length;i++){if(a[i]<b[i])return true;if(a[i]>b[i])return false;}return false;};
  let leaves=0,nodes=0;
  function rec(lo,left) {
   nodes++;
   for(let j=0;j<h;j++)if(sums[j]+left*maxTail[lo][j]<=2*d)return;
   if(left===0) {
    leaves++;
    let g=d;for(const x of a)g=gcd(g,x);if(g!==1)return;
    for(const u of units) {
     const b=a.map(x=>{let r=(u*x)%d;return Math.min(r,d-r)}).sort((x,y)=>x-y);
     if(lexless(b,a))return;
    }
    candidates.push([...a]);return;
   }
   for(let x=lo;x<=h;x++){
    a.push(x);for(let j=0;j<h;j++)sums[j]+=weights[x][j];
    rec(x,left-1);
    for(let j=0;j<h;j++)sums[j]-=weights[x][j];a.pop();
   }
  }
  rec(0,n);byOrder[d]={candidates,nodes,leaves};
 }
 output[n]=byOrder;
}
console.log(JSON.stringify(output));

