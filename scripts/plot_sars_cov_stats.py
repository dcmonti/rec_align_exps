import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def compute_tools_performances(results):
    means_df = pd.DataFrame(columns=["tool", "mean_edit_score", "var_edit_score"])
    tools = results["tool"].unique()
    for tool in tools:
        results_tool = results[results["tool"] == tool]
        mean_edit_score = round(results_tool["edit_score"].mean(),2 )
        var_edit_score = round(results_tool["edit_score"].std(),2 )
        
        means_df = means_df._append({"tool": tool, "mean_edit_score": mean_edit_score, "var_edit_score": var_edit_score}, ignore_index=True)
    means_df.to_csv("output/sars-cov-2/mean_edit_scores.csv", index=False)
    return means_df

def remap_tools(results):
    # use more describing names for tools
    results["tool"] = results["tool"].replace("ga", "GraphAligner")
    results["tool"] = results["tool"].replace("ra", "RecAlign")
    results["tool"] = results["tool"].replace("mc", "Minichain")

def check_correct_path(data):
    # iter on the records
    for index, row in data.iterrows():
        read_name = row["read"]
        haplo, _, _ = read_name.split("_")
        read_type = row["type"]
        paths = row["paths"]
        
        original_haplo = ""
        # get ref path
        if read_type == "ori":
            ref_path = "output/sars-cov-2/sim"
            with open(f"{ref_path}/filtered_sd_{haplo[1:]}_0001.ref", "r") as f:
                for line in f:
                    if line.startswith(">"):
                        original_haplo = line.strip().split(" ")[0][1:]
                        break
            if pd.isna(paths):
                data.loc[index, "correct_path"] = False
            elif original_haplo in paths:
                data.loc[index, "correct_path"] = True
            else:
                data.loc[index, "correct_path"] = False

# Load the data
data = pd.read_csv("output/sars-cov-2/switches_edit_scores.csv")
check_correct_path(data)
remap_tools(data)


# plot edit scores for each tool
sns.boxplot(hue=data["tool"], x="tool", y="edit_score", data=data, showfliers=False)
plt.title("Edit scores for each tool")

plt.savefig("output/sars-cov-2/edit_scores.png")
plt.show()

#plot time for each tool
sns.boxplot(hue=data["tool"], x="tool", y="time", data=data, showfliers=False)
plt.title("Time for each tool")
plt.yscale("log")
plt.savefig("output/sars-cov-2/time.png")
plt.show()  

#plot memory for each tool
sns.boxplot(hue=data["tool"], x="tool", y="memory", data=data, showfliers=False)        
plt.title("Memory for each tool")   
plt.yscale("log")
plt.savefig("output/sars-cov-2/memory.png")
plt.show()

# save mean edit scores
mean_edit_scores = compute_tools_performances(data)
print(mean_edit_scores)

data["error_rate"] = data["edit_score"]/data["read_length"]
# assign category every 0.005 increase in error rate
data["error_rate"] = data["error_rate"].apply(lambda x: round(x*200)/200)
# plot time for each tool, grouped by error rate
sns.boxplot(hue=data["tool"], x="error_rate", y="time", data=data, showfliers=False)
plt.show()

#table with correct path for each tool (omly for original reads)
data = data[data["type"] == "ori"]
correct_paths = data[data["correct_path"] == True]
correct_paths = correct_paths[["tool", "correct_path"]]
correct_paths = correct_paths.groupby("tool").count()
correct_paths = correct_paths.rename(columns={"correct_path": "correct_paths"})
correct_paths["total_reads"] = data.groupby("tool").count()["read"]
#correct_paths["correct_paths"] = correct_paths["correct_paths"]/correct_paths["total_reads"]
print(correct_paths)
