package quant;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class CalDif {
	
	public static void main(String[] arg) throws IOException{
		CalDif test=new CalDif();
		String path="e:/wolfei/20141106-quantArea/allwolfei.out2";
		test.process(path, path.replace("out2", "out3"));
	}
	
	public void process(String path,String out) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine();
		int allrep=7; int[] eachrep={2,3,2};
		while (line!=null){
			String [] array=line.replace("\t\t", "\t-\t-").split("\t");
			String pid=array[0].replace("-", "");
			String pepnumst="";
			for (int i=0;i<allrep;i++){
				pepnumst+=array[i+1].replace("-", "")+"\t";
			}
			String AveErr="";
			ArrayList<Double>list=new ArrayList<>(); //ave, stdev, ave, stdev, ave, stdev
			int col=allrep;
			for (int i=0;i<eachrep.length;i++){
				col+=eachrep[i]+3;
				System.out.println(array.length+" "+col);
				String avest=array[col-2].replace("-", "");
				String stdevst=array[col-1].replace("-", "");
				String pcerror=array[col].replace("-", "");
				//changed line start
				String avelog=avest;
				if (avelog.length()>0)
					avelog=Math.log10(Double.parseDouble(avest))+"";
				AveErr+=avelog+"\t"+pcerror+"\t";
				//changed line end
				double ave=0.0,stdev=0.0;
				if (stdevst.length()>0){
					ave=Double.parseDouble(avest);
					stdev=Double.parseDouble(stdevst);
				}
				list.add(ave);list.add(stdev);
			}
			bw.write(pid+"\t"+pepnumst+AveErr+calDif(list)+"\n");
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}

	private String calDif(ArrayList<Double>list){
		String returns="";
		int sample=list.size()/2;
		for (int i=0;i<list.size()-2;i+=2){
			if (list.get(i)>0 && list.get(i+2)>0){
				double difference=list.get(i)-list.get(i+2);
				double sumstdev=list.get(i+1)+list.get(i+3);
				double dif=(double) Math.round(difference/sumstdev*10)/10;
				//changed
				returns+=Math.log10(Math.abs(difference))+"\t"+dif+"\t";
//				returns+=difference+"\t"+dif+"\t";
			}
			else
				returns+="\t\t";
		}
		if (list.size()>4){
			if (list.get(list.size()-2)>0 && list.get(0)>0){
			double difference=list.get(list.size()-2)-list.get(0);
			double sumstdev=list.get(list.size()-1)+list.get(1);
			double dif=(double) Math.round(difference/sumstdev*10)/10;
			//changed
			returns+=Math.log10(Math.abs(difference))+"\t"+dif+"\t";
			}else
				returns+="\t\t";
		}
		returns=returns.substring(0,returns.length()-1);
		return returns;
	}
}
