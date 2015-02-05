package quant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeMap;

public class aveProt {
	public static void main(String[] arg) throws IOException{
		aveProt test=new aveProt();
//		String path="e:/wolfei/20141106-quantAera/ecoli/allecoli.out";
		String path="e:/wolfei/20141106-quantArea/allwolfei.out";
		test.process(path, path+"2");
	}

	public void process(String path,String out) throws IOException{
		int samplenum=7; double maxarea=0;double maxnorm=0;
		ArrayList<Double> sumlist=new ArrayList<>();
		ArrayList<Double> normlist=new ArrayList<Double>();
		HashSet<String> protID=new HashSet<String>();
		protID.add("114565715");protID.add("114567371");protID.add("114565705");
		for (int i=0;i<samplenum;i++){
			double norm=getNormProt(path,i*4+3,0,protID); //String path,int areacol,int pepcol,int protcol
			double sum=calSum(path,i*4+3,i*4+1,0); //String path,int areacol,int pepcol,int protcol
			if (norm>maxnorm)
				maxnorm=norm;
//			if (sum>maxarea)
//				maxarea=sum;
			sumlist.add(sum);
			normlist.add(norm);
		}
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.replace("\t\t", "\t-\t-").split("\t");
			String pid=array[0].replace("-", "");
			if (!pid.contains("P")&&!pid.contains("Q")&&!pid.contains(";")&&!pid.contains("spike")){
				ArrayList<Double> list=new ArrayList<>();
				ArrayList<String> peplist=new ArrayList<>();
				for (int i=0;i<samplenum;i++){
					String pepst=array[i*4+1].replace("-", "");
					double area=0.0;
					peplist.add(array[i*4+1].replace("-", ""));
					if (pepst.length()>0){
						int pepnu=Integer.parseInt(pepst);
						if (pepnu>2){
							String areast=array[i*4+3].replace("-", "");
							if (areast.length()>0)
								area=Double.parseDouble(array[i*4+3]);
						}
					}			
					list.add(area);
				}
				int[] bpos={0,3};
				int[] cpos={1,4,5};
				int[] mpos={2,6};
				ArrayList<int[]> listpos=new ArrayList<>();
				listpos.add(bpos); listpos.add(cpos);listpos.add(mpos);
				String pepstr="",areastr="";
				for (int i=0;i<listpos.size();i++){
					ArrayList<Double> blist=new ArrayList<>();
					int[] pos=listpos.get(i);
					for (int j=0;j<pos.length;j++){
//						System.out.println(j);
						blist.add(list.get(pos[j]));
						blist.add(normlist.get(pos[j]));
//						blist.add(sumlist.get(pos[j]));

						pepstr+=peplist.get(pos[j])+"\t";
					}
					areastr+=normalizearea(blist,maxnorm);
//					areastr+=normalizearea(blist,maxarea);
				}
				bw.write(array[0].replace("-", "")+"\t"+pepstr+areastr+"\n");	
			}
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}

	private String normalizearea(ArrayList<Double> list,double norm){
		String returns="";
		double sum=0;
		ArrayList<Double>newlist=new ArrayList<>(); //not count 0
		ArrayList<Double> alllist=new ArrayList<>(); //count 0
		for (int i=0;i<list.size()-1;i+=2){ //i: area, i+1: sum area
			double area=list.get(i);		
			double areasum=list.get(i+1);
			double newarea=area*norm/areasum; //change to area of norm prot in highest sample/area of norm prot in this sample
			if (area>0){
				newlist.add(newarea);
				sum+=newarea;
			}
			alllist.add(newarea);
		}

		double ave=sum/newlist.size();
		double sumdifference=0;
		for (int i=0;i<newlist.size();i++){
			double difference=newlist.get(i)-ave;
			sumdifference+=difference*difference;
//			sumdifference=difference*difference;
		}
		for (int i=0;i<alllist.size();i++){
			if (alllist.get(i)>0)
				returns+=alllist.get(i)+"\t";
			else
				returns+="\t";
		}
		if (newlist.size()>1){
			double stdev=(Math.sqrt(sumdifference/(newlist.size()-1)));
//			int pcerror=(int)Math.round(Math.sqrt(sumdifference/(newlist.size()-1))/sum*100);
			int pcerror=(int)Math.round(stdev/ave*100);

			returns+=ave+"\t"+stdev+"\t"+pcerror+"\t";
		}else if (newlist.size()==1)
			returns+=ave+"\t\t\t";
		else
			returns+="\t\t\t";
		return returns;
	}

	//at least 2 pep
	private double calSum(String path,int areacol,int pepcol,int protcol) throws IOException{
		double sum=0;
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.replace("\t", "\t-").split("\t");
			String pep=array[pepcol].replace("-", "");
			String area=array[areacol].replace("-", "");
			String prot=array[protcol].replace("-", "");
			if (!prot.contains("P")&&!prot.contains("Q")&&!prot.contains(";")){
				if (pep.length()>0){
					int pepnum=Integer.parseInt(pep);
					if (pepnum>2){
						if (area.length()>0){
							double areanum=Double.parseDouble(area);
							sum+=areanum;
						}
					}
				}
			}
			line=br.readLine();
		}
		br.close();
		return sum;
	}
	
	private double getNormProt(String path,int areacol,int protcol,HashSet<String> protID) throws IOException{
		double norm=0;
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		double sumarea=0;
		while (line!=null){
			String[] array=line.replace("\t", "\t-").split("\t");
			String area=array[areacol].replace("-", "");
			String prot=array[protcol].replace("-", "");
			if (protID.contains(prot)){
				sumarea+=Double.parseDouble(area);
			}
			line=br.readLine();
		}
		br.close();
		norm=sumarea/protID.size();
		return norm;	
	}

}
