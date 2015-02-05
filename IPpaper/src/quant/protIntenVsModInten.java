package quant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

public class protIntenVsModInten {

	public static void main(String[] arg) throws IOException{
		String path="e:/wolfei/20141010/BCmCBC-SGmod.ptt";
		String modtype="Acetyl,";
		String out="e:/wolfei/20141021/protvsmod-b1-"+modtype+"-modnum.out";
		modnum(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+"-modnum.out";
		modnum(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+"-modnum.out";
		modnum(path,out,21,40,2,49,18,modtype);

		modtype="Butyryl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+"-modnum.out";
		modnum(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+"-modnum.out";
		modnum(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+"-modnum.out";
		modnum(path,out,21,40,2,49,18,modtype);

		modtype="Crotonyl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+"-modnum.out";
		modnum(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+"-modnum.out";
		modnum(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+"-modnum.out";
		modnum(path,out,21,40,2,49,18,modtype);

		modtype="3OHButyryl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+"-modnum.out";
		modnum(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+"-modnum.out";
		modnum(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+"-modnum.out";
		modnum(path,out,21,40,2,49,18,modtype);

		modtype="acetoacetyl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+"-modnum.out";
		modnum(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+"-modnum.out";
		modnum(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+"-modnum.out";
		modnum(path,out,21,40,2,49,18,modtype);
	}
	
	public static void main0(String[] arg) throws IOException{
		String path="e:/wolfei/20141010/BCmCBC-SGmod.ptt";
		String modtype="Acetyl,";
		String out="e:/wolfei/20141021/protvsmod-b1-"+modtype+".out";
		read(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+".out";
		read(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+".out";
		read(path,out,21,40,2,49,18,modtype);

		modtype="Butyryl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+".out";
		read(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+".out";
		read(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+".out";
		read(path,out,21,40,2,49,18,modtype);

		modtype="Crotonyl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+".out";
		read(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+".out";
		read(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+".out";
		read(path,out,21,40,2,49,18,modtype);

		modtype="3OHButyryl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+".out";
		read(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+".out";
		read(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+".out";
		read(path,out,21,40,2,49,18,modtype);

		modtype="acetoacetyl,";
		out="e:/wolfei/20141021/protvsmod-b1-"+modtype+".out";
		read(path,out,13,40,2,41,10,modtype);

		out="e:/wolfei/20141021/protvsmod-c1-"+modtype+".out";
		read(path,out,17,40,2,45,14,modtype);

		out="e:/wolfei/20141021/protvsmod-mc1-"+modtype+".out";
		read(path,out,21,40,2,49,18,modtype);
	}
	
	public static void modnum(String path,String out, int protlogcol,int modcol,int PIDcol, int modpepareacol,int pepnumcol,String modtype) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine(); line=br.readLine();
		TreeMap<String, String> prot=new TreeMap<>(); //pid, log prot area
		TreeMap<String, Integer> modpep=new TreeMap<>(); //pid, area of mod pep
		while (line!=null){
			line=line.replace("\t\t", "\t-\t-");
			String[] array=line.split("\t");
			int len=(array.length);
			String protlog=array[protlogcol];
			String PID=array[PIDcol];
			int pepnum=0;
			String pepnumst=array[pepnumcol].replace("-", "");
			if (!pepnumst.equals(""))
				pepnum=Integer.parseInt(pepnumst);
			if (pepnum>2){
				prot.put(PID, protlog);
				if (modcol<len){
					String mod=array[modcol].replace("-", "");			
					double modpeparea=0.0;
					String modpepareast=array[modpepareacol].replace("-", "");
					if (!modpepareast.equals(""))
						modpeparea=Double.parseDouble(modpepareast);
					if (mod.equals(modtype)){
						if (modpeparea>0){
							if (modpep.containsKey(PID)){
								modpep.put(PID, 1+modpep.get(PID));
							}else
								modpep.put(PID, 1);
						}
					}
				}
			}
			line=br.readLine();
		}
		bw.write("PID\tlog(prot area)\tsumlog(modpeparea)"+modtype+"\n");
		for (Map.Entry<String, String> entry:prot.entrySet()){
			String pid=entry.getKey();
			String logprotarea=entry.getValue();
			if (modpep.containsKey(pid)){
				bw.write(pid+"\t"+logprotarea+"\t"+modpep.get(pid)+"\n");
			}else
				bw.write(pid+"\t"+logprotarea+"\t\t"+"\n");
		}
		br.close();bw.flush();bw.close();
	}

	public static void read(String path,String out, int protlogcol,int modcol,int PIDcol, int modpepareacol,int pepnumcol,String modtype) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine(); line=br.readLine();
		TreeMap<String, String> prot=new TreeMap<>(); //pid, log prot area
		TreeMap<String, Double> modpep=new TreeMap<>(); //pid, area of mod pep
		while (line!=null){
			line=line.replace("\t\t", "\t-\t-");
			String[] array=line.split("\t");
			int len=(array.length);
			String protlog=array[protlogcol];
			String PID=array[PIDcol];
			int pepnum=0;
			String pepnumst=array[pepnumcol].replace("-", "");
			if (!pepnumst.equals(""))
				pepnum=Integer.parseInt(pepnumst);
			if (pepnum>2){
				prot.put(PID, protlog);
				if (modcol<len){
					String mod=array[modcol].replace("-", "");			
					double modpeparea=0.0;
					String modpepareast=array[modpepareacol].replace("-", "");
					if (!modpepareast.equals(""))
						modpeparea=Double.parseDouble(modpepareast);
					if (mod.equals(modtype)){
						if (modpeparea>0){
							if (modpep.containsKey(PID)){
								modpep.put(PID, modpeparea+modpep.get(PID));
							}else
								modpep.put(PID, modpeparea);
						}
					}
				}
			}
			line=br.readLine();
		}
		bw.write("PID\tlog(prot area)\tsumlog(modpeparea)"+modtype+"\n");
		for (Map.Entry<String, String> entry:prot.entrySet()){
			String pid=entry.getKey();
			String logprotarea=entry.getValue();
			if (modpep.containsKey(pid)){
				double logmodarea=Math.log10(modpep.get(pid));
				bw.write(pid+"\t"+logprotarea+"\t"+logmodarea+"\n");
			}else
				bw.write(pid+"\t"+logprotarea+"\t\t"+"\n");
		}
		br.close();bw.flush();bw.close();
	}

}
