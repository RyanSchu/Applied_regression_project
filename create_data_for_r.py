import sys
import pandas as pd 

def read_in_dataframe(df):
    return_df = df.transpose()
    return_df.columns = return_df.iloc[0]
    return_df = return_df[1:]
    return return_df

if __name__ == '__main__':
    target = pd.read_csv(sys.argv[1], sep='\t')
    target = read_in_dataframe(target)
    
    predictors = pd.read_csv(sys.argv[2], sep='\t')
    predictors = read_in_dataframe(predictors)

    # Create the final Product
    joined_df = predictors.merge(target, left_index=True, right_index=True)
    print("Joined dataframes together in pandas")
    print("Outputting to file " + sys.argv[3])
    joined_df.to_csv(sys.argv[3], sep=",", index=False)
