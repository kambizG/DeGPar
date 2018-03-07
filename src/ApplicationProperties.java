import java.io.File;
import java.io.FileInputStream;
import java.util.Properties;

public class ApplicationProperties {

	// Constants
	public String DIRECTORY;
	public String SECTION_MARKER;
	public String SUB_DIRECTORY;

	// Files
	public File GRAPH;

	// Parameters
	public double MINIMUM_WEIGHT_THRESHOLD;
	public int NUMBER_OF_PARTITIONS;
	public double EXCHANGE_LIMIT;
	public double EXPANSION_LIMIT;
	public double TERMINATION_THRESHOLD;
	public int COUNT_STABILITY;
	public String COLORING;
	public boolean PRINT_LOG;
	public int NUM_REPETITION;

	public ApplicationProperties(String appProp, String in_directory) {

		try {
			Properties props = new Properties();
			FileInputStream in = new FileInputStream(appProp);
			props.load(in);
			in.close();
			MINIMUM_WEIGHT_THRESHOLD = Double.parseDouble(props
					.getProperty("MINIMUM_WEIGHT_THRESHOLD"));
			TERMINATION_THRESHOLD = Double.parseDouble(props
					.getProperty("TERMINATION_THRESHOLD"));
			EXCHANGE_LIMIT = Double.parseDouble(props
					.getProperty("EXCHANGE_LIMIT"));
			EXPANSION_LIMIT = Double.parseDouble(props
					.getProperty("EXPANSION_LIMIT"));
			PRINT_LOG = Boolean.parseBoolean(props.getProperty("PRINT_LOG"));
			NUMBER_OF_PARTITIONS = Integer.parseInt(props
					.getProperty("NUMBER_OF_PARTITIONS"));
			NUM_REPETITION = Integer.parseInt(props
					.getProperty("NUM_REPETITION"));
			COUNT_STABILITY = Integer.parseInt(props
					.getProperty("COUNT_STABILITY"));
			SECTION_MARKER = props.getProperty("SECTION_MARKER");
			COLORING = props.getProperty("COLORING");
			DIRECTORY = in_directory;

			GRAPH = new File(DIRECTORY + props.getProperty("GRAPH"));

		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void initiateSubDirectory(int number) {
		SUB_DIRECTORY = DIRECTORY + "trial_" + number + "/";
		//System.currentTimeMillis() + "_"
		//		+ MINIMUM_WEIGHT_THRESHOLD + "_" + NUMBER_OF_PARTITIONS + "_"
		//		+ EXCHANGE_LIMIT + "_" + EXPANSION_LIMIT + "_"
		//		+ TERMINATION_THRESHOLD + "/";
	}

}
