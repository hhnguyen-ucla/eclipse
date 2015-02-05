package quant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

public class modStoi {

	public static void main(String[] arg) throws IOException{
		modStoi test=new modStoi();
		String path="E:/wolfei/20141106-quantArea/mod-013015.all";
		String out=path+"stoi";
//		test.highStoi(path, out);
		String wolfeip="e:/wolfei/20141106-quantArea/20150201prot.txt";
//		test.sumStoi(out, out+"2", wolfeip);
		BufferedWriter bw = new BufferedWriter(new FileWriter("e:/temp.txt"));
		bw.write(""); bw.flush();bw.close();
		test.protStoi(path, path+"prot", wolfeip);
	}
	/*
E:\wolfei\20141106-quantArea\quantArea2\20140928wolfeiB.mod4
E:\wolfei\20141106-quantArea\quantArea2\20141024wolfeiB.mod4
E:\wolfei\20141106-quantArea\quantArea2\20140928wolfeiC.mod4
E:\wolfei\20141106-quantArea\quantArea2\20141024wolfeiC.mod4
E:\wolfei\20141106-quantArea\quantArea2\20141024wolfeiC_2.mod4
E:\wolfei\20141106-quantArea\quantArea2\20140930wolfeimC.mod4
E:\wolfei\20141106-quantArea\quantArea2\20141024wolfeimC.mod4
E:\wolfei\20141106-quantArea\mod-013015.all
	 */	
	//if stoichiometry>10, write 1, or -2, orelse write "". Take merge of .mod4
	public void highStoi(String path,String out) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.replace("\t", "\t-").split("\t");
			int initial=1;
			int groupAcyl=6;
			int groupSubstrate=7;
			int loop=7;
			String stoi="";
			for (int i=0;i<groupSubstrate;i++){
				for (int j=1;j<groupAcyl;j++){
					//col0: title, col6,12,18,24,30: b
					//col31:empty,col37,43,49,55,61
					//col62:empty,col68
					int position=initial-1+(i*5+j)*groupAcyl+i;
//					System.out.println(position);
					String ratiost=array[position].replace("-", "");
					String istoi="";
					if (ratiost.length()>0){
						int ratio=Integer.parseInt(ratiost);
						if (ratio>=10)
							istoi="1";
						else if (ratio<2)
							istoi="2";
						else
							istoi="0";
					}else{
						if (array[position-4].length()>1){
//							System.out.println(array[position-4]);
							istoi="0";
						}
					}
					stoi+="\t"+istoi;
				}			

			}
			bw.write(line+stoi+"\n");
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}
	
	//Y if at least 2 columns highisotopic=1, col (1,2) (3,4,5) (6,7)
	//sum for each growth condition
	//sum for each mod
	//sum for all growth
	//proteins (blast, location, log area, compare)
	public void sumStoi(String path,String out,String wolfeip) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine();
		TreeMap<String, String>tree=readwolfeiProt(wolfeip);
		while (line!=null){
			String[] array=line.replace("\t", "\t@").split("\t");
			String[] oldkey=array[0].split(",");
			String pid=oldkey[0];
			String value=tree.get(pid);
			int start=array.length-7*5;
			String[] highstoi=new String[15];
			int pos=0;
			for (int i=start;i<start+5;i++){
				String[] b=new String[]{array[i].replace("@", ""),array[i+5].replace("@", "")};
				System.out.println(i+" "+(int)(i+5)+" "+determineStoi(b));
				highstoi[pos]=determineStoi(b);
				pos++;
			}
			for (int i=start+10;i<start+15;i++){
				String[] c=new String[]{array[i].replace("@", ""),array[i+5].replace("@", ""),array[i+10].replace("@", "")};
				System.out.println(i+" "+(int)(i+5)+(int)(i+10));
				highstoi[pos]=determineStoi(c);
				pos++;
			}
			for (int i=start+25;i<start+30;i++){
				String[]m=new String[]{array[i].replace("@", ""),array[i+5].replace("@", "")};
				highstoi[pos]=determineStoi(m);
				System.out.println(i+" "+(int)(i+5));

				pos++;
			}
			int[] sumeachgrowth=new int [3];
			int[] loweachgrowth=new int[3];
			int[] loweachmod=new int [5];
			int[] sumeachmod=new int [5];
			for (int i=0;i<3;i++){
				int sum=0;int low=0;
				for (int j=0;j<5;j++){
					String st=highstoi[5*i+j];
					if (st.equals("H"))
						sum++;
					if (st.equals("L"))
						low++;
					
				}
				sumeachgrowth[i]=sum;
				loweachgrowth[i]=low;
			}
			for (int i=0;i<5;i++){
				int sum=0;int low=0;
				for (int j=0;j<3;j++){
					//0,5,10; 1,6,11;
					String st=highstoi[j*5+i];
					if (st.equals("H"))
						sum++;
					if (st.equals("L"))
						low++;
				}
				sumeachmod[i]=sum; loweachmod[i]=low;
			}			
//			int[] sumall=new int[3];//high,low,all
			int high=0,low=0,nd=0;
			for (int i=0;i<highstoi.length;i++){
				String st=highstoi[i];
				if (st.equals("H"))
					high++;
				if (st.equals("L"))
					low++;
				if (st.equals("ND"))
					nd++;
			}
			nd=nd+high+low;
			String highst=high+"",lowst=low+"",all=nd+"";
			if (high==0)
				highst="";
			if (low==0)
				lowst="";
			if (nd==0)
				all="";
			
			bw.write(pid+value+array[0].replace(pid, "").replace(",", "\t")+"\t");
			bw.write(highst+"\t"+lowst+"\t"+all+"\t"); //high,low,all
			for (int i=0;i<sumeachmod.length;i++){
				String write=sumeachmod[i]+"";
				if (sumeachmod[i]==0)
					write="";
				bw.write(write+"\t");
			}
			for (int i=0;i<loweachmod.length;i++){
				String write=loweachmod[i]+"";
				if (loweachmod[i]==0)
					write="";
				bw.write(write+"\t");
			}
			for (int i=0;i<sumeachgrowth.length;i++){
				String write=sumeachgrowth[i]+"";
				if (sumeachgrowth[i]==0)
					write="";
				bw.write(write+"\t");
			}
			for (int i=0;i<loweachgrowth.length;i++){
				String write=loweachgrowth[i]+"";
				if (loweachgrowth[i]==0)
					write="";
				bw.write(write+"\t");
			}			
			for (int i=0;i<highstoi.length;i++){
				bw.write(highstoi[i]+"\t");
			}
			bw.write(line.replace(array[0]+"\t", "")+"\n");
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}
	private String determineStoi(String[] stoi ){
		//H if at least 2 are "1"
		//L if at least 2 are "-2"
		//"" if all are ""
		//ND if not H or L and at least 1 are not ""
		String returns="ND";
		int h=0,l=0,none=0;
		for (int i=0;i<stoi.length;i++){
			if (stoi[i].equals("1"))
				h++;
			else if (stoi[i].equals("2"))
				l++;
			else if (stoi[i].equals(""))
				none++;
		}
		if (h>1)
			returns="H";
		if (l>1)
			returns="L";
		if (none==stoi.length)
			returns="";
		return returns;
		
	}
	
	private TreeMap<String, String>readwolfeiProt(String path) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		TreeMap<String, String> tree=new TreeMap<String, String>();
		String line=br.readLine();
		while (line!=null){
			String[] array=line.split("\t");
			String key=array[0]; //pid
			String value=line.replace(key, "");
			tree.put(key, value);
			line=br.readLine();
		}
		br.close();
		return tree;
	}
	
	
	//protein (blast, location, log area, compare), sum high low each growth, sum high low all growth
	public void protStoi(String path,String out,String wolfeip) throws IOException{
		TreeMap<String, String>tree=readwolfeiProt(wolfeip);
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine();
		TreeMap<String, int[]>mods=new TreeMap<String, int[]>();
		while (line!=null){
			String[] array=line.replace("\t", "\t@").split("\t");
			String[] oldkey=array[0].split(",");
			String pid=oldkey[0];
			int[] combine=(processaline(line));
			if (mods.containsKey(pid)){
				int[] vmods=mods.get(pid);
				for (int i=0;i<vmods.length;i++){
					vmods[i]=vmods[i]+combine[i];
				}
				mods.put(pid, vmods);
			}else
				mods.put(pid,combine);
			line=br.readLine();
		}
		for (Map.Entry<String, int[]> e:mods.entrySet()){
			String pid=e.getKey();
			int[] vmods=e.getValue(); //6slots*5 mod*3 growth
			String value=tree.get(pid);
			tree.remove(pid);
			bw.write(pid+value+"\t");
			for (int i=0;i<vmods.length;i+=6){
//				System.out.println(i);
				int high=vmods[i]+vmods[i+3];
				int low=vmods[i+1]+vmods[i+4];
				int total=0;
				String sum="";
				for (int x=i;x<i+6;x++){
					total+=vmods[x];
					sum+=vmods[x]+":";
				}
				String highst=high+"",lowst=low+"",totalst=total+"";
				if (high==0) highst="";
				if (low==0) lowst="";
				if (total==0) totalst="";
				sum=sum.substring(0,sum.length()-1);
				bw.write(highst+"\t"+lowst+"\t"+totalst+"\t"+sum+"\t");
			}
			bw.write("\n");
		}		
		
		for (Map.Entry<String, String> e:tree.entrySet()){
			bw.write(e.getKey()+e.getValue()+"\t\t\t\t\t\n");
		}
		br.close();bw.flush();bw.close();
	}
	
	private int[] combinealine(int[] count){		
//		System.out.println("count len "+count.length);
		int[] returns=new int[90]; //high,low,notder,high*,low*,notder*: 6*5mods*3growth
		for (int i=0;i<returns.length;i++)
			returns[i]=0;
		String[] highstoi=new String[15];
		for (int j=0;j<5;j++){
			int[] a=new int[6];
			int[] b=new int[6];
			for (int i=0;i<6;i++){
				a[i]=count[i+6*j]; //0->5,6->11,12->17,18->23,24->29
				//last set: j=4, a:i+24:24->29				
				b[i]=count[i+6*j+30];//30->59
				//last set: j=4, b:5+24+30=
			}
			int[] c=determinestoi2(a, b);
			for (int i=0;i<6;i++)
				returns[i+6*j]=c[i]; //0->29
		}
		
		for (int j=0;j<5;j++){
			int[] a=new int[6];
			int[] b=new int[6];
			int[] c=new int [6];
			for (int i=60;i<66;i++){
//				System.out.println(i+6*j);
				a[i-60]=count[i+6*j]; //60->89
				b[i-60]=count[i+6*j+30];//90->119
				c[i-60]=count[i+6*j+60];//120->149: 65+24+60=149
			}
			int[] d=determinestoi3(a, b,c);
			for (int i=30;i<36;i++)
				returns[i+6*j]=d[i-30]; //30->59
		}
		for (int j=0;j<5;j++){
			int[] a=new int[6];
			int[] b=new int[6];
			for (int i=150;i<156;i++){
				a[i-150]=count[i+6*j]; //0->6,7->13,14->20,21->27,28->34
				b[i-150]=count[i+6*j+30];//35->41,42->48,
			}
			int[] c=determinestoi2(a, b);
			for (int i=60;i<66;i++)
				returns[i+6*j]=c[i-60];
		}
		return returns;
	}
	
	private int[] determinestoi2(int[] a,int[]b){
		int[] returns=new int[6]; //h,l,n,h*,l*,n*
		for (int i=0;i<returns.length;i++)
			returns[i]=0;
		for (int i=0;i<2;i++){
			if (a[i]+b[i]>1)
				returns[i]=1; //0,1,3,4 pairs
			if (a[i+3]+b[i+3]>1)
				returns[i+3]=1;
		}
		for (int i=0;i<2;i++){
			if (a[i]+b[i+3]>1 || a[i+3]+b[i]>1) //0-3, 1-4 pairs
				returns[i]=1;
		}
		int sum=0;
		for (int i=0;i<returns.length;i++){
			sum+=returns[i];
		}
		int orisum=0;
		for (int i=0;i<a.length;i++){
			orisum+=a[i]+b[i];
		}
		if (sum<1 && orisum>0){
			boolean sharepep=true;
			for (int i=0;i<3;i++){
				if (a[i]>0 || b[i]>0)
					sharepep=false;
			}
			if (sharepep==false)
				returns[2]=1;
			else
				returns[5]=1;
		}
		return returns;		
	}
	
	private int[] determinestoi3(int[] a,int[]b,int[]c){
		int[] returns=new int[6]; //h,l,n,h*,l*,n*
		for (int i=0;i<returns.length;i++)
			returns[i]=0;
		for (int i=0;i<2;i++){
			if (a[i]+b[i]+c[i]+a[3]+b[3]+c[3]>1){
				if (a[i]>0||b[i]>0||c[i]>0)
					returns[i]=1;
				else
					returns[i+3]=1;
			}
		}
		int sum=0;
		for (int i=0;i<returns.length;i++){
			sum+=returns[i];
		}
		int orisum=0;
		for (int i=0;i<a.length;i++){
			orisum+=a[i]+b[i]+c[i];
		}
		if (sum<1&&orisum>0){
			boolean sharepep=true;
			for (int i=0;i<3;i++){
				if (a[i]>0 || b[i]>0||c[i]>0)
					sharepep=false;
			}
			if (sharepep==false)
				returns[2]=1;
			else
				returns[5]=1;
		}
		return returns;		
	}
	
	private int[] processaline(String line) throws IOException{
		int[] count=new int[210]; //high,low,notder,high*,low*,notder*: 6*5mod*7samples
		for (int i=0;i<count.length;i++)
			count[i]=0;
		String[] array=line.replace("\t", "\t@").split("\t");
		
		int initial=1;
		int groupAcyl=6;
		int groupSubstrate=7;
		String stoi="";
		int countpos=0;
		for (int i=0;i<groupSubstrate;i++){
			for (int j=1;j<groupAcyl;j++){
				int position=initial-1+(i*5+j)*groupAcyl+i;
				String ratiost=array[position].replace("@", "");
				String istoi="";
				if (ratiost.length()>0){
					int ratio=Integer.parseInt(ratiost);
					if (ratio>=10)
						istoi="1";
					else if (ratio<2)
						istoi="2";
					else
						istoi="0";
				}else{
					if (array[position-4].length()>1){
						istoi="0";
					}
				}
				int shift=0;
				String scorelen=array[position-4];
				if (scorelen.contains("*")){
					boolean toshift=true;
					String[] a=scorelen.split(":");
					for (int k=0;k<a.length;k++){
						if (!a[k].contains("*"))
							toshift=true;
					}
					if (toshift==true)
						shift=3;
				}
				if (istoi.equals("1"))
					count[countpos+shift+0]=1;
				if (istoi.equals("2"))
					count[countpos+shift+1]=1;
				if (istoi.equals("0"))
					count[countpos+shift+2]=1;				
				countpos+=6;				
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("e:/temp.txt",true));
		int[] combine=combinealine(count);
		bw.write(line+"\t");
		for (int i=0;i<count.length;i++){
			if (count[i]<1)
				bw.write("\t");
			else
				bw.write(count[i]+"\t");
		}
		for (int i=0;i<combine.length;i++){
			if (combine[i]<1)
				bw.write("\t");
			else
				bw.write(combine[i]+"\t");
		}
		bw.write("\n");
		bw.flush();bw.close();
		return combine;
		
	}
}
