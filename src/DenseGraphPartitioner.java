import java.io.File;

public class DenseGraphPartitioner {

	static Logger logger;

	public static void main(String[] args) {

		if(args.length == 0){
			System.out.println("ERR: Not Enough Parameters...");
			System.exit(0);
		}
		String directory = args[0];
		//String directory = "experiments/exp_09/";
		try {
			// Load Properties
			// -----------------------------------------------------------
			ApplicationProperties ap = new ApplicationProperties("appProperties.properties", directory);
			// Partitioning-----------------------------------------------------------------
			for (int i = 0; i < ap.NUM_REPETITION; i++) {
				ap.initiateSubFolder(i);
				File f = new File(ap.SUB_FOLDER);
				if (!f.exists())
					if (f.mkdir());
				logger = new Logger(ap.SUB_FOLDER + "report");
				GraphStruct g = new GraphStruct(ap, logger);
				long start = System.currentTimeMillis();
				g.start();
				logger.logEntry("################");
				logger.logEntry("#LOG: Execution terminated in: "
						+ (System.currentTimeMillis() - start) / 1000.0 + "s");
				g.writeEdgeList(ap.SUB_FOLDER + "edgeList.csv", 0.0);
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}