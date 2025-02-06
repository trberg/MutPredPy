
import os
import pandas as pd
import scipy.io as sio


class Catalog:

    def __init__(self, job_dir, storage_location, dry_run):
        
        self.job_dir = job_dir
        self.storage_location = storage_location
        self.dry_run = dry_run


    def create_index(self):
        pass


    def check_storage_location_setup(self):

        
        if os.path.exists(self.storage_location) and os.path.exists(f"{self.storage_location}/catalog"):
            pass
        else:
            os.makedirs(self.storage_location)
            os.makedirs(f"{self.storage_location}/catalog")


        if os.path.exists(f"{self.storage_location}/index.json"):
            pass
        else:
            self.create_index()



    def load_features(self, job):

        feature_list = "/Users/bergqt01/Research/MutPredPy/MutPredPy/resources/mp2_features.VPcats.csv"
        features = pd.read_csv(feature_list)

        scores = pd.read_csv(f"{self.job_dir}/{job}/output.txt")

        cur_file = 1
        starting_index = 0
        
        all_features = []

        while os.path.exists(f"{self.job_dir}/{job}/output.txt.feats_{cur_file}.mat"):
            
            mat_file = sio.loadmat(f"{self.job_dir}/{job}/output.txt.feats_{cur_file}.mat")

            indexed_scores = scores.iloc[starting_index: starting_index + mat_file["feats"].shape[0]]
            
            multiindex = pd.MultiIndex.from_frame(indexed_scores[["ID","Substitution"]], names=["ID","Substitution"])
            
            cur_file += 1
            starting_index += mat_file["feats"].shape[0]

            features_data = pd.DataFrame(data=mat_file["feats"], columns=list(features["Feature name"]), index=multiindex)
            features_data = features_data.reset_index()

            all_features.append(features_data)
        
        all_features = pd.concat(all_features)
        
        return all_features
        



    def load_mutpred_file(self, job, file_type):

        property_list = "/Users/bergqt01/Research/MutPredPy/MutPredPy/resources/sandbox_features_ordering.csv"
        properties = pd.read_csv(property_list, header=None, names=["Feature name"])

        feature_list = "/Users/bergqt01/Research/MutPredPy/MutPredPy/resources/mp2_features.VPcats.csv"
        features = pd.read_csv(feature_list)

        different_names = [properties, features]

        names = {len(L):L for L in different_names}

        scores = pd.read_csv(f"{self.job_dir}/{job}/output.txt")

        cur_file = 1
        starting_index = 0
        
        all_props = []
        
        while os.path.exists(f"{self.job_dir}/{job}/output.txt.{file_type}_{cur_file}.mat"):
            
            mat_file = sio.loadmat(f"{self.job_dir}/{job}/output.txt.{file_type}_{cur_file}.mat")

            indexed_scores = scores.iloc[starting_index: starting_index + mat_file[f"{file_type}"].shape[0]]
            multiindex = pd.MultiIndex.from_frame(indexed_scores[["ID","Substitution"]], names=["ID","Substitution"])

            prop_data = pd.DataFrame(data=mat_file[f"{file_type}"], columns=list(names[mat_file[f"{file_type}"].shape[1]]["Feature name"]), index=multiindex)

            all_props.append(prop_data)

            cur_file += 1
            starting_index += 1


        all_props = pd.concat(all_props)

        return all_props
    


    def info(self, job):

        print (os.listdir(f"{self.job_dir}/{job}"))

        mat_file = sio.loadmat(f"{self.job_dir}/{job}/output.txt.propX_pu_1.mat")
        
        #print (mat_file)
        print (mat_file["propX_pu"])
        print (mat_file["propX_pu"].shape)

        print (mat_file["propX_pu"][0])
        
        #scores = pd.read_csv(f"{self.job_dir}/{job}/output.txt")
        #print (scores)




if __name__=="__main__":

    cat = Catalog("/Users/bergqt01/Research/Scratch/Testing/jobs", "/Users/bergqt01/Research/Scratch/Testing/mutpred2_catalog", dry_run=True)
    
    property_list = "/Users/bergqt01/Research/MutPredPy/MutPredPy/resources/sandbox_features_ordering.csv"
    properties = pd.read_csv(property_list, header=None, names=["Feature name"])
    
    properties = list(properties["Feature name"])

    features = cat.load_mutpred_file("job_203", "feats")
    #print ([col for col in features.columns if col.split("_")[0] in properties])
    #print (len([col for col in features.columns if col.split("_")[0] in properties]))
    #exit()
    #print (features[[col for col in features.columns if col.split("_")[0] in properties["Feature name"]]])
    #exit()
    print (features[[col for col in features.columns if "VSL2B_disorder" in col]])
    

    #posts = cat.load_posteriors("job_99")
    #print (posts)

    #prop_scores = cat.load_mutpred_file("job_203", "prop_scores_pu")
    #print (prop_scores[[col for col in prop_scores.columns if "Proteolytic" in col]])

    cat.info("job_203")

