% Group the particles of same type into subgroups by checking if they are connected to each other
function Group = CONNECTIVITY2(Category, neighboridx)
Group = {};
while ~isempty(Category)
    
    % create a new list using the first particle in remaining list of particle indices
    size_old = 0;
    size_new = 1;
    List = Category(1);
    
    % add current list particles' neighbors to the list if they have same type. Opt out the loop
    % when the list does not expand
    while size_old ~= size_new && length(Group) < 5
        new_entry_idx = List(size_old+1:end);
        size_old = size_new;
        for new_pt = new_entry_idx
            for j = neighboridx{new_pt}
                if ismember(j,Category) && ~ismember(j,List)
                    List = cat(2,List,j);
                end
            end
        end
        size_new = size(List,2);
    end
    count = zeros(1,size(List,2));
    for i = 1:size(List,2)
        for j = neighboridx{List(i)}
            if ismember(j,List)
                count(i) = count(i)+1;
            end
        end
    end
    if any(count>4) 
        List(count<=3) = [];
    end
    temp = CONNECTIVITY(List, neighboridx);
    for i = 1:length(temp)
        for j = temp{i}
            Category(Category == j) = [];
        end
        Group =cat(2,Group,sort(temp{i}));
    end
end

for i = 1:size(Group,2)
    group_size(i) = length(Group{i});
end
[~,idx] = sort(group_size,'descend');
Group = Group(idx);
end