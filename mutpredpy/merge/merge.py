import pandas as pd
import os


class Merge:
    def __init__(self, output, job_dir, mechanisms, dry_run):

        self.output = output
        self.job_dir = job_dir
        self.mechanisms = mechanisms
        self.dry_run = dry_run
        self.scored = self.mutpred_scores()

    def write_output(self, data, output_file):
        data.to_csv(output_file, sep="\t", index=False)

    def mutpred_scores(self):

        if self.mechanisms:
            included_columns = [
                "ID",
                "Substitution",
                "MutPred2 score",
                "Molecular mechanisms with Pr >= 0.01 and P < 1.00",
                "Motif information",
                "Remarks",
            ]

        else:
            included_columns = ["ID", "Substitution", "MutPred2 score"]

        print("MERGE")
        scores = pd.concat(
            [
                pd.read_csv(f"{self.job_dir}/{m}/output.txt", usecols=included_columns)
                for m in os.listdir(self.job_dir)
                if os.path.exists(f"{self.job_dir}/{m}/output.txt")
            ]
        )

        scores = scores.drop_duplicates()

        return scores

    def merge(self):

        self.write_output(self.scored, self.output)
